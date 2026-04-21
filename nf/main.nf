nextflow.enable.dsl = 2

def mustExist = { label, p ->
    if (!p) {
        error "Missing required parameter: ${label}"
    }
    def f = file(p.toString())
    if (!f.exists()) {
        error "Not found: ${label}=${p} (launchDir=${launchDir}, projectDir=${projectDir})"
    }
    f.toString()
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
    sample_col <- "Source Name"

    if (is.na(pop_col) || is.na(rna_col) || !sample_col %in% names(sdrf)) {
      stop("Could not detect required SDRF columns.")
    }

    targets <- unlist(strsplit("${eur_pop_arg}", ",", fixed=TRUE))
    sdrf_eur <- sdrf[get(pop_col) %in% targets]

    mapped <- merge(sdrf_eur, fam, by.x=sample_col, by.y="V2")
    mapped <- mapped[!duplicated(V2)]

    write.table(mapped[, .(V1, get(sample_col))], "eur_keep.txt",
                sep="\\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    map_out <- mapped[, .(RNA_ID=get(rna_col), IID=get(sample_col), FID=V1)]
    fwrite(map_out, "eur_map.tsv", sep="\\t")

    cat("Mapped", nrow(map_out), "European samples.\\n")
    """
}

process GENERATE_PCA {
    publishDir "${params.outdir}/pca", mode: 'copy'

    input:
    path bed
    path bim
    path fam
    path keep_list

    output:
    path "geno_pca.eigenvec", emit: eigenvec
    path "genotype_grm*",     emit: grm

    script:
    """
    ${params.gcta_path} --bfile ${bed.baseName} \
        --keep ${keep_list} \
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

    output:
    path "final.phenotypes.tsv", emit: pheno
    path "final.qcovar", emit: qcovar
    path "gene_index_map.txt", emit: map
    path "filtered_gene_ids.txt", emit: filtered_genes

    script:
    """
    echo "Starting phenotype prep: TPM+count filter (GTEx-style mandatory), TMM on counts, then INT"

    Rscript ${params.prepare_pheno_script} \
      ${tpm} \
      ${counts} \
      ${id_map} \
      ${pca_file} \
      ${fam_file} \
      final \
      ${num_peer} \
      ${params.gtex_tpm_threshold} \
      ${params.gtex_count_threshold} \
      ${params.gtex_sample_frac_threshold} \
      ${params.n_pcs}
    """
}

process ESTIMATE_HERITABILITY {
    tag "${gene_name}"
    errorStrategy 'ignore'
    maxRetries 0

    input:
    tuple val(gene_name), val(idx)
    path pheno
    path qcovar
    path grm_files

    output:
    path "${gene_name}.stats", emit: stats

    script:
    """
    ${params.gcta_path} --grm ${params.grm_prefix} \
        --pheno ${pheno} --mpheno ${idx} --qcovar ${qcovar} \
        --reml --out ${gene_name}_greml --thread-num 1 || echo "GREML_FAIL"

    ${params.gcta_path} --grm ${params.grm_prefix} \
        --pheno ${pheno} --mpheno ${idx} --qcovar ${qcovar} \
        --HEreg --out ${gene_name}_he --thread-num 1 || echo "HE_FAIL"

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

    echo -e "${gene_name}\\t\${status}\\t${idx}\\t\${h2_reml}\\t\${se_reml}\\t\${pval_reml}\\t\${h2_he}\\t\${se_he}" > ${gene_name}.stats
    """
}

process SUMMARIZE_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path all_stats

    output:
    path "final_heritability_summary.tsv"

    script:
    """
    echo -e "Gene\\tStatus\\tIndex\\th2_GREML\\tSE_GREML\\tPval_GREML\\th2_HE\\tSE_HE" > final_heritability_summary.tsv
    cat *.stats >> final_heritability_summary.tsv
    """
}

workflow {
    log.info "launchDir=${launchDir}"
    log.info "projectDir=${projectDir}"
    log.info "workDir=${workflow.workDir}"
    log.info "outdir=${params.outdir}"

    params.tpm_file = mustExist("tpm_file", params.tpm_file)
    params.counts_file = mustExist("counts_file", params.counts_file)
    params.prepare_pheno_script = mustExist("prepare_pheno_script", params.prepare_pheno_script)
    def bed_path = mustExist("plink bed", "${params.plink_prefix}.bed")
    def bim_path = mustExist("plink bim", "${params.plink_prefix}.bim")
    def fam_path = mustExist("plink fam", "${params.plink_prefix}.fam")

    tpm_ch = Channel.value(file(params.tpm_file))
    counts_ch = Channel.value(file(params.counts_file))
    bed_ch = Channel.value(file(bed_path))
    bim_ch = Channel.value(file(bim_path))
    fam_ch = Channel.value(file(fam_path))

    FILTER_EUROPEANS(params.sdrf_url, fam_ch, params.eur_pops)

    GENERATE_PCA(bed_ch, bim_ch, fam_ch, FILTER_EUROPEANS.out.keep_file)

    PREPARE_PHENOTYPES(
        tpm_ch,
        counts_ch,
        FILTER_EUROPEANS.out.map_file,
        GENERATE_PCA.out.eigenvec,
        fam_ch,
        params.peer_nk
    )

    gene_ch = PREPARE_PHENOTYPES.out.map
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.gene_name, row.mpheno_index) }

    ESTIMATE_HERITABILITY(
        gene_ch,
        PREPARE_PHENOTYPES.out.pheno.collect(),
        PREPARE_PHENOTYPES.out.qcovar.collect(),
        GENERATE_PCA.out.grm.collect()
    )

    SUMMARIZE_RESULTS(ESTIMATE_HERITABILITY.out.stats.collect())
}
