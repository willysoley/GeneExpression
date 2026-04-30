import java.security.MessageDigest
import java.time.Instant

nextflow.enable.dsl = 2

def mustExist = { label, p ->
    if (!p) {
        error "Missing required parameter: ${label}"
    }
    def f = file(p.toString())
    if (!f.exists()) {
        error "Not found: ${label}=${p} (launchDir=${launchDir}, projectDir=${projectDir})"
    }
    f.toAbsolutePath().normalize().toString()
}

def requiredPrepareDiagnostics = [
    "Applying mandatory GTEx-style filter",
    "Running TMM normalization on raw counts",
    "Expression source:",
    "Normalization mode:",
    "Phenotype prep completed successfully"
]

def sha256Hex = { pathStr ->
    def md = MessageDigest.getInstance("SHA-256")
    new File(pathStr.toString()).withInputStream { is ->
        byte[] buffer = new byte[8192]
        for (int read = is.read(buffer); read != -1; read = is.read(buffer)) {
            md.update(buffer, 0, read)
        }
    }
    md.digest().collect { String.format("%02x", it) }.join()
}

def validatePreparePhenoScript = { label, scriptPath, requiredMarkers ->
    def resolvedPath = mustExist(label, scriptPath)
    def scriptFile = new File(resolvedPath)
    def scriptText = scriptFile.getText("UTF-8")
    def missingMarkers = requiredMarkers.findAll { marker -> !scriptText.contains(marker) }
    def mtimeUtc = Instant.ofEpochMilli(scriptFile.lastModified()).toString()
    def sha256 = sha256Hex(resolvedPath)

    log.info "${label}=${resolvedPath}"
    log.info "${label}_mtime_utc=${mtimeUtc}"
    log.info "${label}_sha256=${sha256}"

    if (missingMarkers) {
        error "${label} missing required diagnostic marker(s): ${missingMarkers.join(' | ')}"
    }
    log.info "${label} diagnostics preflight OK (${requiredMarkers.size()} markers)"

    [path: resolvedPath, sha256: sha256, mtimeUtc: mtimeUtc]
}

def normalizeMaybeQuotedEmpty = { value ->
    def s = value?.toString()
    if (s == null) return ""
    s = s.trim()
    if (s == '""' || s == "''") return ""
    s
}

process FILTER_EUROPEANS {
    input:
    val sdrf_url
    path fam_file
    val eur_pops

    output:
    path "eur_keep.txt", emit: keep_file
    path "eur_map.tsv",  emit: map_file

    script:
    def eur_pop_arg = eur_pops.join(",")
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
      library(data.table)
    })

    if (!file.exists("sdrf.txt")) download.file("${sdrf_url}", "sdrf.txt")

    sdrf <- fread("sdrf.txt", header=TRUE, sep="\\t", check.names=FALSE)
    names(sdrf) <- make.unique(names(sdrf))
    fam  <- fread("${fam_file}", header=FALSE)

    pop_col <- grep("ancestry category", names(sdrf), ignore.case=TRUE, value=TRUE)[1]
    rna_col <- grep("ENA_RUN", names(sdrf), value=TRUE)[1]

    if (is.na(pop_col) || is.na(rna_col)) {
      stop("Could not detect required SDRF columns: ancestry category and/or ENA_RUN.")
    }

    targets <- unlist(strsplit("${eur_pop_arg}", ",", fixed=TRUE))
    sdrf_eur <- sdrf[get(pop_col) %in% targets]

    fam[, V2 := as.character(V2)]
    candidate_sample_cols <- intersect(
      c("Source Name", "Characteristics[individual]", "Factor Value[individual]"),
      names(sdrf_eur)
    )
    if (length(candidate_sample_cols) == 0L) {
      stop("No candidate sample-ID column found in SDRF.")
    }

    overlap_n <- vapply(
      candidate_sample_cols,
      function(col) sum(as.character(sdrf_eur[[col]]) %in% fam[["V2"]], na.rm=TRUE),
      numeric(1)
    )
    sample_col <- candidate_sample_cols[which.max(overlap_n)]
    if (is.na(sample_col) || overlap_n[which.max(overlap_n)] == 0) {
      stop("Could not map SDRF sample IDs to FAM IID (V2).")
    }

    mapped <- merge(sdrf_eur, fam, by.x=sample_col, by.y="V2")
    mapped <- unique(mapped, by=sample_col)

    write.table(mapped[, .(V1, get(sample_col))], "eur_keep.txt",
                sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    map_out <- mapped[, .(RNA_ID=get(rna_col), IID=get(sample_col), FID=V1)]
    fwrite(map_out, "eur_map.tsv", sep="\\t")

    cat("Mapped", nrow(map_out), "European samples.\\n")
    """
}

process GENERATE_PCA {
    publishDir { "${params.outdir}/pca/${pca_label}" }, mode: 'copy'

    input:
    path bed
    path bim
    path fam
    path keep_list
    val hm3_extract_arg
    val pca_label

    output:
    path "geno_pca.eigenvec", emit: eigenvec
    path "genotype_grm*",     emit: grm

    script:
    """
    ${params.gcta_path} --bfile ${bed.baseName} \
        --keep ${keep_list} \
        ${hm3_extract_arg} \
        --make-grm --out genotype_grm --thread-num ${task.cpus}

    ${params.gcta_path} --grm genotype_grm --pca ${params.n_pcs} --out geno_pca --thread-num ${task.cpus}
    """
}

process PREPARE_PHENOTYPES {
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    path tpm
    path counts
    path id_map
    path pca_file
    path fam_file
    val num_peer
    val peer_max_genes
    val pheno_prefix

    output:
    path "*.phenotypes.tsv", emit: pheno, arity: '1'
    path "*.qcovar", emit: qcovar, arity: '1'
    path "*.gene_index_map.txt", emit: map, arity: '1'
    path "*.filtered_gene_ids.txt", emit: filtered_genes, arity: '1'

    script:
    def runtimeDiagnostics = requiredPrepareDiagnostics
        .collect { marker -> "\"${marker}\"" }
        .join("\n      ")
    def fallbackPrefix = "pheno_${params.use_hm3_no_hla ? 'hm3_no_mhc' : 'all_snps'}_${params.expression_source}_${params.normalization_type}"
    def incomingPhenoPrefix = pheno_prefix?.toString() ?: ""
    """
    set -euo pipefail
    pheno_prefix_resolved="${incomingPhenoPrefix}"
    if [[ -z "\${pheno_prefix_resolved}" ]]; then
      pheno_prefix_resolved="${fallbackPrefix}"
      echo "WARN: Received empty pheno_prefix in PREPARE_PHENOTYPES; using fallback=\${pheno_prefix_resolved}"
    fi

    echo "Starting phenotype prep: TPM+count filter (GTEx-style mandatory), expression_source=${params.expression_source}, normalization=${params.normalization_type}"

    Rscript ${params.prepare_pheno_script} \
      "${tpm}" \
      "${counts}" \
      "${id_map}" \
      "${pca_file}" \
      "${fam_file}" \
      "\${pheno_prefix_resolved}" \
      "${num_peer}" \
      "${params.gtex_tpm_threshold}" \
      "${params.gtex_count_threshold}" \
      "${params.gtex_sample_frac_threshold}" \
      "${params.n_pcs}" \
      "${params.expression_source}" \
      "${params.normalization_type}" \
      "${peer_max_genes}" \
      2>&1 | tee prepare_phenotypes.log

    required_diagnostics=(
      ${runtimeDiagnostics}
    )
    for marker in "\${required_diagnostics[@]}"; do
      if ! grep -Fq "\${marker}" prepare_phenotypes.log; then
        echo "ERROR: Missing required prepare_phenotypes diagnostic line: \${marker}" >&2
        echo "Checked script: ${params.prepare_pheno_script}" >&2
        exit 1
      fi
    done
    if ! grep -Fq "Normalization mode: ${params.normalization_type}" prepare_phenotypes.log; then
      echo "ERROR: Normalization mode did not propagate to prepare_phenotypes.R" >&2
      exit 1
    fi
    if ! grep -Fq "Expression source: ${params.expression_source}" prepare_phenotypes.log; then
      echo "ERROR: Expression source did not propagate to prepare_phenotypes.R" >&2
      exit 1
    fi

    expected_outputs=(
      "\${pheno_prefix_resolved}.phenotypes.tsv"
      "\${pheno_prefix_resolved}.qcovar"
      "\${pheno_prefix_resolved}.gene_index_map.txt"
      "\${pheno_prefix_resolved}.filtered_gene_ids.txt"
    )
    for f in "\${expected_outputs[@]}"; do
      if [[ ! -s "\${f}" ]]; then
        echo "ERROR: Expected phenotype output missing or empty: \${f}" >&2
        ls -lah . >&2
        exit 1
      fi
    done
    """
}

process ESTIMATE_HERITABILITY {
    tag "${gene_name}"
    errorStrategy { task.exitStatus in [104, 130, 134, 137, 139, 140, 143, 247, 271] ? 'retry' : 'ignore' }
    maxRetries 2
    publishDir "${params.outdir}/diagnostics/heritability", mode: 'copy'

    input:
    tuple val(gene_name), val(idx)
    path pheno
    path qcovar
    path grm_files

    output:
    path "${gene_name}.stats", emit: stats
    path "${gene_name}.diagnostics.tsv", emit: diagnostics
    path "${gene_name}_*.gcta.log", emit: gcta_logs, optional: true

    script:
    // GRM files are staged via `grm_files`; GCTA reads them by this fixed basename.
    def grmPrefix = "genotype_grm"
    """
    set -uo pipefail

    stats_file="${gene_name}.stats"
    diag_file="${gene_name}.diagnostics.tsv"
    greml_log="${gene_name}_greml.gcta.log"
    he_log="${gene_name}_he.gcta.log"

    # Seed outputs early so interrupted tasks still leave breadcrumbs for debugging.
    echo -e "${gene_name}\\tFAIL\\t${idx}\\tNA\\tNA\\tNA\\tNA\\tNA" > "\${stats_file}"
    {
      echo -e "key\\tvalue"
      echo -e "gene\\t${gene_name}"
      echo -e "index\\t${idx}"
      echo -e "start_utc\\t\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      echo -e "workdir\\t\${PWD}"
      echo -e "grm_prefix\\t${grmPrefix}"
      echo -e "pheno_file\\t${pheno}"
      echo -e "qcovar_file\\t${qcovar}"
    } > "\${diag_file}"

    on_term() {
      local sig="\$1"
      {
        echo -e "terminated_signal\\t\${sig}"
        echo -e "terminated_utc\\t\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      } >> "\${diag_file}"
      echo -e "${gene_name}\\tFAIL\\t${idx}\\tNA\\tNA\\tNA\\tNA\\tNA" > "\${stats_file}"
      exit 0
    }
    trap 'on_term TERM' TERM
    trap 'on_term INT' INT

    ${params.gcta_path} --grm ${grmPrefix} \
        --pheno ${pheno} --mpheno ${idx} --qcovar ${qcovar} \
        --reml --out ${gene_name}_greml --thread-num 1 \
        > "\${greml_log}" 2>&1
    greml_rc=\$?

    ${params.gcta_path} --grm ${grmPrefix} \
        --pheno ${pheno} --mpheno ${idx} --qcovar ${qcovar} \
        --HEreg --out ${gene_name}_he --thread-num 1 \
        > "\${he_log}" 2>&1
    he_rc=\$?

    if [ -f "${gene_name}_greml.hsq" ]; then
        h2_reml=\$(grep "V(G)/Vp" ${gene_name}_greml.hsq | awk '{print \$2}')
        se_reml=\$(grep "V(G)/Vp" ${gene_name}_greml.hsq | awk '{print \$3}')
        pval_reml=\$(grep "Pval" ${gene_name}_greml.hsq | awk '{print \$2}')
        status="PASS"
    else
        h2_reml="NA"; se_reml="NA"; pval_reml="NA"; status="FAIL"
    fi

    if [ -f "${gene_name}_he.heri" ]; then
        h2_he=\$(grep "V(G)/Vp" ${gene_name}_he.heri | awk '{print \$2}')
        se_he=\$(grep "V(G)/Vp" ${gene_name}_he.heri | awk '{print \$3}')
    else
        h2_he="NA"; se_he="NA"
    fi

    if [ "\${greml_rc}" -ne 0 ] || [ "\${he_rc}" -ne 0 ]; then
      status="FAIL"
    fi

    {
      echo -e "greml_rc\\t\${greml_rc}"
      echo -e "he_rc\\t\${he_rc}"
      echo -e "greml_hsq_exists\\t\$([ -s ${gene_name}_greml.hsq ] && echo 1 || echo 0)"
      echo -e "he_heri_exists\\t\$([ -s ${gene_name}_he.heri ] && echo 1 || echo 0)"
      echo -e "end_utc\\t\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } >> "\${diag_file}"

    echo -e "${gene_name}\\t\${status}\\t${idx}\\t\${h2_reml}\\t\${se_reml}\\t\${pval_reml}\\t\${h2_he}\\t\${se_he}" > "\${stats_file}"
    """
}

process SUMMARIZE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path all_stats
    val summary_file

    output:
    path "*.tsv", emit: summary, arity: '1'

    script:
    """
    summary_file_resolved="${summary_file}"
    if [[ -z "\${summary_file_resolved}" ]]; then
      summary_file_resolved="final_heritability_summary.tsv"
      echo "WARN: summary_file was empty; using fallback=\${summary_file_resolved}" >&2
    fi
    echo -e "Gene\\tStatus\\tIndex\\th2_GREML\\tSE_GREML\\tPval_GREML\\th2_HE\\tSE_HE" > "\${summary_file_resolved}"
    cat *.stats >> "\${summary_file_resolved}"
    test -s "\${summary_file_resolved}"
    """
}

process COMPRESS_PHENOTYPE {
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    path pheno_tsv
    path summary_file

    output:
    path "*.phenotypes.tsv.gz", emit: pheno_gz

    script:
    """
    set -euo pipefail
    test -s "${summary_file}"
    gzip -c "${pheno_tsv}" > "${pheno_tsv.getName()}.gz"
    """
}

workflow {
    log.info "launchDir=${launchDir}"
    log.info "projectDir=${projectDir}"
    log.info "workDir=${workflow.workDir}"
    log.info "outdir=${params.outdir}"

    params.tpm_file = mustExist("tpm_file", params.tpm_file)
    params.counts_file = mustExist("counts_file", params.counts_file)
    params.use_hm3_no_hla = params.use_hm3_no_hla in [true, "true", "TRUE", "1", 1]
    params.expression_source = (params.expression_source ?: "tmm").toString().trim().toLowerCase()
    def allowedExpressionSource = ["tmm", "tpm"]
    if (!allowedExpressionSource.contains(params.expression_source)) {
        error "Invalid expression_source='${params.expression_source}'. Allowed: ${allowedExpressionSource.join(', ')}"
    }
    params.normalization_type = (params.normalization_type ?: "irnt").toString().trim().toLowerCase()
    if (params.normalization_type == "ukb_irnt") {
        log.info "normalization_type=ukb_irnt is deprecated; using irnt."
        params.normalization_type = "irnt"
    }
    if (params.normalization_type == "tmm_only") {
        log.info "normalization_type=tmm_only is deprecated; using raw."
        params.normalization_type = "raw"
    }
    def allowedNormalization = ["irnt", "inverse_normal", "raw"]
    if (!allowedNormalization.contains(params.normalization_type)) {
        error "Invalid normalization_type='${params.normalization_type}'. Allowed: ${allowedNormalization.join(', ')}"
    }
    log.info "Phenotype expression_source=${params.expression_source}"
    log.info "Phenotype normalization_type=${params.normalization_type}"
    log.info "Phenotype peer_max_genes=${params.peer_max_genes}"
    def snpSetLabel = params.use_hm3_no_hla ? "hm3_no_mhc" : "all_snps"
    def runLabel = "${snpSetLabel}_${params.expression_source}_${params.normalization_type}"
    def phenotypePrefix = normalizeMaybeQuotedEmpty(params.phenotype_prefix)
    if (!phenotypePrefix) {
        params.phenotype_prefix = "pheno_${runLabel}"
    } else {
        params.phenotype_prefix = phenotypePrefix
    }
    log.info "Phenotype prefix=${params.phenotype_prefix}"
    def summaryFilename = normalizeMaybeQuotedEmpty(params.summary_filename)
    if (!summaryFilename) {
        params.summary_filename = "final_heritability_summary_${runLabel}.tsv"
    } else {
        params.summary_filename = summaryFilename
    }
    log.info "Summary filename=${params.summary_filename}"

    def canonicalPrepCandidates = [
        "${projectDir}/nf/bin/prepare_phenotypes.R",
        "${projectDir}/bin/prepare_phenotypes.R",
        "${launchDir}/nf/bin/prepare_phenotypes.R"
    ].unique()
    def canonicalPrepScript = canonicalPrepCandidates.find { file(it).exists() }
    if (!canonicalPrepScript) {
        error "Could not locate canonical nf/bin/prepare_phenotypes.R. Checked: ${canonicalPrepCandidates.join(', ')}"
    }
    def canonicalPrepMeta = validatePreparePhenoScript(
        "canonical_prepare_pheno_script",
        canonicalPrepScript,
        requiredPrepareDiagnostics
    )

    def prepDefaults = [
        "${projectDir}/nf/bin/prepare_phenotypes.R",
        "${projectDir}/bin/prepare_phenotypes.R",
        "${launchDir}/nf/bin/prepare_phenotypes.R"
    ].unique()
    def prepScript = params.prepare_pheno_script?.toString()
    if (!prepScript || !file(prepScript).exists() || prepDefaults.contains(prepScript)) {
        prepScript = prepDefaults.find { file(it).exists() } ?: prepScript
    }
    def selectedPrepMeta = validatePreparePhenoScript(
        "prepare_pheno_script",
        prepScript,
        requiredPrepareDiagnostics
    )
    params.prepare_pheno_script = selectedPrepMeta.path
    if (selectedPrepMeta.sha256 != canonicalPrepMeta.sha256) {
        error "prepare_pheno_script SHA-256 (${selectedPrepMeta.sha256}) differs from canonical nf/bin copy (${canonicalPrepMeta.sha256}). Refusing to run stale or unexpected script."
    }
    def bed_path = mustExist("plink bed", "${params.plink_prefix}.bed")
    def bim_path = mustExist("plink bim", "${params.plink_prefix}.bim")
    def fam_path = mustExist("plink fam", "${params.plink_prefix}.fam")
    def hm3_extract_arg = ""
    if (params.use_hm3_no_hla) {
        def hm3_path = mustExist("hm3_no_hla_snplist", params.hm3_no_hla_snplist)
        hm3_extract_arg = "--extract ${hm3_path}"
        log.info "SNP mode: HM3 no-HLA enabled via ${hm3_path}"
    } else {
        log.info "SNP mode: all variants from plink_prefix (HM3 no-HLA disabled)"
    }

    def tpm_ch = Channel.value(file(params.tpm_file))
    def counts_ch = Channel.value(file(params.counts_file))
    def bed_ch = Channel.value(file(bed_path))
    def bim_ch = Channel.value(file(bim_path))
    def fam_ch = Channel.value(file(fam_path))

    FILTER_EUROPEANS(params.sdrf_url, fam_ch, params.eur_pops)

    GENERATE_PCA(bed_ch, bim_ch, fam_ch, FILTER_EUROPEANS.out.keep_file, hm3_extract_arg, runLabel)

    def reusePhenoDir = params.reuse_pheno_dir?.toString()?.trim()
    def pheno_ch
    def qcovar_ch
    def map_ch

    if (reusePhenoDir) {
        def reuseDirPath = mustExist("reuse_pheno_dir", reusePhenoDir)
        def prefPheno = file("${reuseDirPath}/${params.phenotype_prefix}.phenotypes.tsv")
        def prefQcovar = file("${reuseDirPath}/${params.phenotype_prefix}.qcovar")
        def prefMap = file("${reuseDirPath}/${params.phenotype_prefix}.gene_index_map.txt")

        def reusedPheno
        def reusedQcovar
        def reusedMap

        if (prefPheno.exists() && prefQcovar.exists() && prefMap.exists()) {
            reusedPheno = mustExist("reused phenotype file", prefPheno.toString())
            reusedQcovar = mustExist("reused qcovar file", prefQcovar.toString())
            reusedMap = mustExist("reused gene index map", prefMap.toString())
            log.info "Phenotype mode: reusing combo-specific files (${params.phenotype_prefix}) from ${reuseDirPath}"
        } else {
            def legacyPheno = mustExist("reused legacy phenotype file", "${reuseDirPath}/final.phenotypes.tsv")
            def legacyQcovar = mustExist("reused legacy qcovar file", "${reuseDirPath}/final.qcovar")
            def legacyMap = mustExist("reused legacy gene index map", "${reuseDirPath}/gene_index_map.txt")
            reusedPheno = legacyPheno
            reusedQcovar = legacyQcovar
            reusedMap = legacyMap
            log.info "Phenotype mode: combo-specific files not found; using legacy final.* files from ${reuseDirPath}"
        }
        pheno_ch = Channel.value(file(reusedPheno))
        qcovar_ch = Channel.value(file(reusedQcovar))
        map_ch = Channel.value(file(reusedMap))
    } else {
        PREPARE_PHENOTYPES(
            tpm_ch,
            counts_ch,
            FILTER_EUROPEANS.out.map_file,
            GENERATE_PCA.out.eigenvec,
            fam_ch,
            params.peer_nk,
            params.peer_max_genes,
            params.phenotype_prefix
        )
        pheno_ch = PREPARE_PHENOTYPES.out.pheno
        qcovar_ch = PREPARE_PHENOTYPES.out.qcovar
        map_ch = PREPARE_PHENOTYPES.out.map
    }

    def gene_ch = map_ch
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.gene_name, row.mpheno_index) }

    ESTIMATE_HERITABILITY(
        gene_ch,
        pheno_ch,
        qcovar_ch,
        GENERATE_PCA.out.grm.collect()
    )

    def summaryFileCh = Channel.value(params.summary_filename)
    def allStatsCh = ESTIMATE_HERITABILITY.out.stats.collect()
    SUMMARIZE_RESULTS(allStatsCh, summaryFileCh)
    COMPRESS_PHENOTYPE(pheno_ch, SUMMARIZE_RESULTS.out.summary)
}
