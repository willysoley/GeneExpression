#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(grid)
})

input_path <- "docs/greml_workflow_knowledge_transfer.md"
output_path <- "docs/greml_workflow_knowledge_transfer.pdf"

page_width <- 8.5
page_height <- 11
margin_left <- 0.65
margin_right <- 0.65
margin_top <- 0.55
margin_bottom <- 0.60
body_width_chars <- 104L
code_width_chars <- 104L

parse_markdown <- function(lines) {
  blocks <- list()
  para <- character(0)
  i <- 1L

  add_block <- function(type, text = character(0), level = NA_integer_) {
    blocks[[length(blocks) + 1L]] <<- list(type = type, text = text, level = level)
  }

  flush_para <- function() {
    if (length(para) > 0L) {
      add_block("para", paste(para, collapse = " "))
      para <<- character(0)
    }
  }

  while (i <= length(lines)) {
    line <- lines[[i]]

    if (grepl("^```", line)) {
      flush_para()
      i <- i + 1L
      code <- character(0)
      while (i <= length(lines) && !grepl("^```", lines[[i]])) {
        code <- c(code, lines[[i]])
        i <- i + 1L
      }
      add_block("code", code)
      i <- i + 1L
      next
    }

    if (identical(line, "")) {
      flush_para()
      add_block("space", "")
      i <- i + 1L
      next
    }

    if (grepl("^# ", line)) {
      flush_para()
      add_block("heading", sub("^# ", "", line), level = 1L)
      i <- i + 1L
      next
    }

    if (grepl("^## ", line)) {
      flush_para()
      add_block("heading", sub("^## ", "", line), level = 2L)
      i <- i + 1L
      next
    }

    if (grepl("^### ", line)) {
      flush_para()
      add_block("heading", sub("^### ", "", line), level = 3L)
      i <- i + 1L
      next
    }

    if (grepl("^- ", line)) {
      flush_para()
      add_block("bullet", sub("^- ", "", line))
      i <- i + 1L
      next
    }

    if (grepl("^[0-9]+\\. ", line)) {
      flush_para()
      prefix <- sub("^([0-9]+\\.).*$", "\\1", line)
      text <- sub("^[0-9]+\\. ", "", line)
      add_block("number", paste(prefix, text, sep = "\t"))
      i <- i + 1L
      next
    }

    para <- c(para, line)
    i <- i + 1L
  }

  flush_para()
  blocks
}

wrap_with_prefix <- function(text, first_prefix, next_prefix, width) {
  usable_width <- width - max(nchar(first_prefix), nchar(next_prefix))
  wrapped <- strwrap(text, width = usable_width)

  if (length(wrapped) == 0L) {
    return(first_prefix)
  }

  c(
    paste0(first_prefix, wrapped[[1]]),
    if (length(wrapped) > 1L) paste0(next_prefix, wrapped[-1]) else character(0)
  )
}

hard_wrap <- function(text, width) {
  if (!nzchar(text) || nchar(text) <= width) {
    return(text)
  }

  starts <- seq(1L, nchar(text), by = width)
  substring(text, starts, pmin(starts + width - 1L, nchar(text)))
}

markdown_block_to_lines <- function(block) {
  if (block$type == "heading") {
    width <- if (block$level == 1L) 58L else if (block$level == 2L) 78L else 92L
    return(strwrap(block$text, width = width))
  }

  if (block$type == "para") {
    return(strwrap(block$text, width = body_width_chars))
  }

  if (block$type == "bullet") {
    return(wrap_with_prefix(block$text, first_prefix = "* ", next_prefix = "  ", width = body_width_chars))
  }

  if (block$type == "number") {
    parts <- strsplit(block$text, "\t", fixed = TRUE)[[1]]
    prefix <- paste0(parts[[1]], " ")
    text <- parts[[2]]
    return(wrap_with_prefix(text, first_prefix = prefix, next_prefix = strrep(" ", nchar(prefix)), width = body_width_chars))
  }

  if (block$type == "code") {
    return(unlist(lapply(block$text, hard_wrap, width = code_width_chars), use.names = FALSE))
  }

  character(0)
}

block_style <- function(block) {
  if (block$type == "heading" && block$level == 1L) {
    return(list(fontfamily = "Helvetica", fontface = "bold", fontsize = 15, lineheight = 1.10, pre = 0.00, post = 0.16))
  }

  if (block$type == "heading" && block$level == 2L) {
    return(list(fontfamily = "Helvetica", fontface = "bold", fontsize = 12.5, lineheight = 1.08, pre = 0.15, post = 0.10))
  }

  if (block$type == "heading" && block$level == 3L) {
    return(list(fontfamily = "Helvetica", fontface = "bold", fontsize = 10.8, lineheight = 1.05, pre = 0.10, post = 0.06))
  }

  if (block$type == "code") {
    return(list(fontfamily = "Courier", fontface = "plain", fontsize = 7.4, lineheight = 1.00, pre = 0.04, post = 0.08))
  }

  if (block$type %in% c("bullet", "number")) {
    return(list(fontfamily = "Helvetica", fontface = "plain", fontsize = 9.0, lineheight = 1.08, pre = 0.01, post = 0.02))
  }

  list(fontfamily = "Helvetica", fontface = "plain", fontsize = 9.2, lineheight = 1.10, pre = 0.02, post = 0.07)
}

line_height_inches <- function(fontsize, lineheight) {
  fontsize * lineheight / 72
}

render_pdf <- function(blocks, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  page_num <- 0L
  cursor <- margin_top
  page_open <- FALSE

  draw_footer <- function() {
    grid.text(
      sprintf("Page %d", page_num),
      x = unit(page_width - margin_right, "inches"),
      y = unit(0.28, "inches"),
      just = c("right", "bottom"),
      gp = gpar(fontfamily = "Helvetica", fontsize = 7.5, col = "grey35")
    )
  }

  start_page <- function() {
    page_num <<- page_num + 1L
    cursor <<- margin_top
    page_open <<- TRUE
    grid.newpage()
  }

  finish_page <- function() {
    if (page_open) {
      draw_footer()
      page_open <<- FALSE
    }
  }

  ensure_space <- function(needed) {
    if (!page_open) {
      start_page()
    }

    if ((cursor + needed) > (page_height - margin_bottom)) {
      finish_page()
      start_page()
    }
  }

  draw_rule <- function(y, col = "#30576d") {
    grid.lines(
      x = unit(c(margin_left, page_width - margin_right), "inches"),
      y = unit(rep(page_height - y, 2L), "inches"),
      gp = gpar(col = col, lwd = 0.8)
    )
  }

  pdf(path, width = page_width, height = page_height, family = "Helvetica")
  on.exit(dev.off(), add = TRUE)

  for (block in blocks) {
    if (block$type == "space") {
      ensure_space(0.08)
      cursor <- cursor + 0.08
      next
    }

    lines <- markdown_block_to_lines(block)
    style <- block_style(block)
    lh <- line_height_inches(style$fontsize, style$lineheight)
    needed <- style$pre + max(length(lines), 1L) * lh + style$post

    if (block$type == "heading") {
      needed <- needed + if (block$level == 2L) 0.05 else 0
    }

    ensure_space(needed)
    cursor <- cursor + style$pre

    if (block$type == "code" && length(lines) > 0L) {
      rect_height <- length(lines) * lh + 0.10
      grid.rect(
        x = unit(margin_left - 0.05, "inches"),
        y = unit(page_height - cursor - rect_height / 2 + 0.02, "inches"),
        width = unit(page_width - margin_left - margin_right + 0.10, "inches"),
        height = unit(rect_height, "inches"),
        just = c("left", "center"),
        gp = gpar(fill = "#f5f7f8", col = "#d9e0e3", lwd = 0.4)
      )
      cursor <- cursor + 0.05
    }

    for (line in lines) {
      grid.text(
        line,
        x = unit(margin_left, "inches"),
        y = unit(page_height - cursor, "inches"),
        just = c("left", "top"),
        gp = gpar(
          fontfamily = style$fontfamily,
          fontface = style$fontface,
          fontsize = style$fontsize,
          lineheight = style$lineheight,
          col = if (block$type == "code") "#24333a" else "#111111"
        )
      )
      cursor <- cursor + lh
    }

    if (block$type == "heading" && block$level == 2L) {
      cursor <- cursor + 0.01
      draw_rule(cursor)
      cursor <- cursor + 0.04
    }

    cursor <- cursor + style$post
  }

  finish_page()
  invisible(page_num)
}

if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path)
}

raw_lines <- readLines(input_path, warn = FALSE)
blocks <- parse_markdown(raw_lines)
page_count <- render_pdf(blocks, output_path)

cat("Wrote PDF to:", normalizePath(output_path), "\n")
cat("Pages:", page_count, "\n")
