
# detach("package:magrittr", unload = TRUE)
library("ggplot2")
library("cowplot")
library("dplyr")

user_path <- Sys.getenv("user_path", unset = Sys.getenv("USER_PATH", unset = "/path/to/user"))

library("gg.gap")
library("reshape2")
library("ggbreak")
library("xlsx")
library("ggpmisc")
# install.packages('magrittr')
# library("magrittr");
# install.packages("dplyr")
library("dplyr")
library("openxlsx")
library("ggforce") # devtools::install_github("thomasp85/ggforce", ref = '4008a2e') 
source(file.path(user_path, "lrft2/scripts", "model.R"))

# color_shape_matrix
caller_dm6 <- c("MELT","TEMP2","TLDR_pass", "TLDR_all" , "TELR", "TrEMOLO",  "LOCATE_pass","LOCATE_all" )
caller_hg38 <- c("MELT","TEMP2",    "TLDR_pass", "TLDR_all","xTea_pass","xTea_all",  "PALMER" ,  "MEHunter", "GraffiTE", "LOCATE_pass", "LOCATE_all" )
caller_dm6_dis <- c("MELT", "TEMP2",  "TLDR_pass", "TELR", "TrEMOLO", "TLDR_pass-LOCATE12-pass", "TLDR_pass-spe", "LOCATE12-pass-TLDR_pass", "LOCATE12-pass-spe", "Olp_3-LOCATE12", "Olp_3-TELR", "Olp_3-TLDR", "TLDR_all", "LOCATE12-all", "LOCATE12-pass")
caller_hg38_dis <- c("TLDR_pass", "LOCATE_pass", "LOCATE_all", "xTea_pass", "xTea_all", "TLDR_pass-LOCATE.new.v2-pass","TLDR_pass-spe","LOCATE.new.v2-pass-TLDR_pass","LOCATE.new.v2-pass-spe")


callers <- c("TEMP2",  "MELT",      "TLDR_pass", "TLDR_all",  "xTea_pass", "xTea_all", "TELR",    "TrEMOLO",    "PALMER",  "MEHunter", "GraffiTE", "LOCATE_all", "LOCATE_pass" ) # , "LOCATE_pass-TLDR_pass", "LOCATE_pass-spe", "TLDR_pass-LOCATE_pass", "TLDR_pass-spe"
# colors <- c("#06DA93",  "#00aa00", "#FFCD00",   "#f6b26b" ,  "#D642CA", "#AEAFAE",  "#008000",   "#459f5e",   "#FF7903", "#2A9D8F",  "#FF0000",     "#CC0000") # ,    "#ff3535",               "#e06666",          "#e69138",              "#ff9900"
colors <- c("#FFC107", "#FF6F00",   "#2196F3",   "#2196F3" ,  "#4CAF50", "#4CAF50",   "#BDBDBD", "#8D6E63",    "#3F51B5",   "#7E57C2",  "#f1b6da",  "red",     "red") # ,    "#ff3535",               "#e06666",          "#e69138",              "#ff9900"
# colors <- c("#FF7043", "#FFC107",   "#2196F3",   "#2196F3" ,  "#4CAF50", "#4CAF50",   "#26A69A", "#26C6DA",    "#red",   "#7E57C2",    "red",     "red") # ,    "#ff3535",               "#e06666",          "#e69138",              "#ff9900"

line_types <- c("solid", "solid",  "solid",      "dashed",  "solid",   "dashed",     "solid",    "solid",   "solid",   "solid",  "solid",   "dashed",      "solid"    )
shapes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
line_size <- c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4)
color_shape_matrix <- data.frame(caller = callers, color=colors, lines=line_types,  shape=shapes, size=line_size)
rownames(color_shape_matrix) <- callers


# read data
genome <- "hg38"
result_root <- file.path(user_path, "2022_long_reads", "result")
data_path <- file.path(result_root, "giab", "intersection")
perfomance_dis <- read.table(paste0(data_path, "/intersection_revision/all.merge.dis.txt"), header = T)
# perfomance_dis <- perfomance_dis[perfomance_dis$frequency>0.2,]
head(perfomance_dis)




perfomance_dis.LOCATE_TLDR_xTea <- perfomance_dis[perfomance_dis$LOCATE.div!="-" & perfomance_dis$TLDR.div!="-" & perfomance_dis$xTea.div!="-", , ]
perfomance.LOCATE_TLDR_xTea.te_len <- melt(perfomance_dis.LOCATE_TLDR_xTea, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE.te_len","TLDR.te_len" ,"xTea.te_len"  ), variable.name = "caller", value.name = "te_len")
perfomance.LOCATE_TLDR_xTea.te_len$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_xTea.te_len$caller)
perfomance.LOCATE_TLDR_xTea.div <- melt(perfomance_dis.LOCATE_TLDR_xTea, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE.div","TLDR.div","xTea.div" ), variable.name = "caller", value.name = "div")
perfomance.LOCATE_TLDR_xTea.div$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_xTea.div$caller)
perfomance.LOCATE_TLDR_xTea.len_div <- merge(perfomance.LOCATE_TLDR_xTea.te_len, perfomance.LOCATE_TLDR_xTea.div, by = c("insID","sample", "depth", "caller", "TE"))
perfomance.LOCATE_TLDR_xTea.len_div <- perfomance.LOCATE_TLDR_xTea.len_div[perfomance.LOCATE_TLDR_xTea.len_div$div != "-", ]
perfomance.LOCATE_TLDR_xTea.len_div$div <- as.numeric(perfomance.LOCATE_TLDR_xTea.len_div$div)
perfomance.LOCATE_TLDR_xTea.len_div$te_len <- as.numeric(perfomance.LOCATE_TLDR_xTea.len_div$te_len)
head(perfomance.LOCATE_TLDR_xTea.len_div)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(xml2)
})

output_dir <- file.path(user_path, "2022_long_reads", "figure_revision")

caller_colors <- c(
  "LOCATE" = "red",
  "TLDR" = "#2196F3",
  "xTea" = "#4CAF50"
)

caller_order <- c("LOCATE", "TLDR", "xTea")
te_order <- c("ALU", "LINE1", "SVA")

default_theme <- theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_line(size = 0.15, colour = "grey88"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 0.25, colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "black"),
    axis.text = element_text(colour = "black", size = 7),
    axis.title = element_text(colour = "black", size = 8),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.25),
    strip.text = element_text(colour = "black", size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.35, "cm"),
    legend.key.height = unit(0.25, "cm"),
    plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
    plot.margin = margin(3, 3, 3, 3, "mm")
  )

plot_theme <- function() {
  if (exists("mytemp_locate", inherits = TRUE)) {
    get("mytemp_locate", inherits = TRUE)
  } else {
    default_theme
  }
}

excel_col_to_int <- function(x) {
  chars <- utf8ToInt(x)
  Reduce(function(total, char) total * 26 + char - utf8ToInt("A") + 1, chars, init = 0)
}

read_xlsx_first_sheet_xml <- function(path) {
  tmp_dir <- tempfile("xlsx_xml_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  utils::unzip(
    path,
    files = c("xl/sharedStrings.xml", "xl/worksheets/sheet1.xml"),
    exdir = tmp_dir
  )

  ns <- c(x = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")

  shared_doc <- xml2::read_xml(file.path(tmp_dir, "xl/sharedStrings.xml"))
  shared_strings <- vapply(xml2::xml_find_all(shared_doc, ".//x:si", ns), function(si) {
    paste0(xml2::xml_text(xml2::xml_find_all(si, ".//x:t", ns)), collapse = "")
  }, character(1))

  sheet_doc <- xml2::read_xml(file.path(tmp_dir, "xl/worksheets/sheet1.xml"))
  rows <- xml2::xml_find_all(sheet_doc, ".//x:sheetData/x:row", ns)

  row_values <- lapply(rows, function(row) {
    cells <- xml2::xml_find_all(row, "./x:c", ns)
    values <- rep(NA_character_, 7)

    for (cell in cells) {
      ref <- xml2::xml_attr(cell, "r")
      col <- excel_col_to_int(gsub("[0-9]", "", ref))
      if (col > length(values)) next

      raw_value <- xml2::xml_text(xml2::xml_find_first(cell, "./x:v", ns))
      if (is.na(raw_value) || raw_value == "") next

      values[col] <- if (identical(xml2::xml_attr(cell, "t"), "s")) {
        shared_strings[as.integer(raw_value) + 1]
      } else {
        raw_value
      }
    }

    values
  })

  header <- row_values[[1]]
  data <- as.data.frame(do.call(rbind, row_values[-1]), stringsAsFactors = FALSE)
  names(data) <- header
  data
}

summarise_div <- function(div_df, group_cols) {
  div_df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_div = mean(div, na.rm = TRUE),
      sd_div = stats::sd(div, na.rm = TRUE),
      se_div = sd_div / sqrt(n),
      .groups = "drop"
    )
}

prepare_div_data <- function(div_df) {
  required_cols <- c("caller", "TE", "div")
  missing_cols <- setdiff(required_cols, names(div_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  div_df %>%
    dplyr::mutate(
      caller = factor(as.character(caller), levels = caller_order),
      TE = factor(as.character(TE), levels = te_order),
      div = as.numeric(div)
    ) %>%
    dplyr::filter(!is.na(caller), !is.na(TE), !is.na(div))
}

plot_div_overall <- function(div_df) {
  plot_data <- div_df %>% summarise_div("caller")

  ggplot(plot_data, aes(x = caller, y = mean_div, fill = caller)) +
    geom_col(width = 0.65, colour = "black", size = 0.2) +
    geom_errorbar(
      aes(ymin = pmax(mean_div - se_div, 0), ymax = mean_div + se_div),
      width = 0.18,
      size = 0.2
    ) +
    scale_fill_manual(values = caller_colors) +
    scale_y_continuous(
      expand = c(0.02, 0),
      labels = function(x) sprintf("%.2f", x)
    ) +
    labs(x = NULL, y = "Mean divergence", title = "GIAB") +
    theme(legend.position = "none") +
    plot_theme()
}

plot_div_by_te <- function(div_df) {
  plot_data <- div_df %>% summarise_div(c("TE", "caller"))

  ggplot(plot_data, aes(x = caller, y = mean_div, fill = caller)) +
    geom_col(width = 0.65, colour = "black", size = 0.2) +
    geom_errorbar(
      aes(ymin = pmax(mean_div - se_div, 0), ymax = mean_div + se_div),
      width = 0.18,
      size = 0.2
    ) +
    scale_fill_manual(values = caller_colors) +
    scale_y_continuous(
      expand = c(0.02, 0),
      labels = function(x) sprintf("%.2f", x)
    ) +
    labs(x = NULL, y = "Mean divergence", title = "GIAB by TE") +
    facet_wrap(~TE, nrow = 1, scales = "free_y") +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    plot_theme()
}

plot_giab_div <- function(div_df, output_prefix = "giab.sequence_divergence") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  div_data <- prepare_div_data(div_df)
  overall_summary <- summarise_div(div_data, "caller")
  by_te_summary <- summarise_div(div_data, c("TE", "caller"))

  write.csv(
    overall_summary,
    file.path(output_dir, paste0(output_prefix, ".overall.summary.csv")),
    row.names = FALSE
  )
  write.csv(
    by_te_summary,
    file.path(output_dir, paste0(output_prefix, ".by_TE.summary.csv")),
    row.names = FALSE
  )

  p_overall <- plot_div_overall(div_data)
  p_by_te <- plot_div_by_te(div_data)

  ggsave(
    file.path(output_dir, paste0(output_prefix, ".overall.pdf")),
    p_overall,
    width = 2.2,
    height = 2.2,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0(output_prefix, ".overall.png")),
    p_overall,
    width = 2.2,
    height = 2.2,
    units = "in",
    dpi = 300
  )
  ggsave(
    file.path(output_dir, paste0(output_prefix, ".by_TE.pdf")),
    p_by_te,
    width = 5.0,
    height = 2.4,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    file.path(output_dir, paste0(output_prefix, ".by_TE.png")),
    p_by_te,
    width = 5.0,
    height = 2.4,
    units = "in",
    dpi = 300
  )

  invisible(list(
    overall_plot = p_overall,
    by_te_plot = p_by_te,
    overall_summary = overall_summary,
    by_te_summary = by_te_summary
  ))
}

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else getwd()
}

input_file <- file.path(dirname(script_path), "giab.div.xlsx")
is_rscript_run <- !interactive() && length(grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)) > 0

if (is_rscript_run) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }
  div_df <- read_xlsx_first_sheet_xml(input_file)
  plot_giab_div(div_df)
}

# Usage in R:
# source("plot_giab_div_bar.R")
# div_df <- read_xlsx_first_sheet_xml("giab.div.xlsx")
# result <- plot_giab_div(div_df)
#
# p1 <- result$overall_plot
# p2 <- result$by_te_plot
div_df <- read_xlsx_first_sheet_xml("giab.div.xlsx")
result <- plot_giab_div(div_df)

p1 <- result$overall_plot
p2 <- result$by_te_plot

pic <- cowplot::plot_grid(p1, p2, nrow = 1, ncol=2)
# pic.1 <- cowplot::plot_grid(pic, pic_all[[1]][[2]],nrow = 2, ncol=1, rel_heights=c(1, 0.4),scale = c(1,1))
# Save Excel file.

ggsave(pic, filename = file.path(output_dir, "giab.sequence_divergence.pdf"), width = 6, height = 2.4)
