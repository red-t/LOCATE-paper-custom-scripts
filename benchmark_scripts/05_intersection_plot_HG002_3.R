
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


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

# args <- commandArgs(trailingOnly = FALSE)
# script_arg <- "--file="
# script_path <- sub(script_arg, "", args[grep(script_arg, args)])
# script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path)) else getwd()


sample <- "HG002_7"

result_root <- file.path(user_path, "2022_long_reads", "result")
input_file <- file.path(result_root, sample, "intersection_revision", "pp_intersection", "data", "performance.txt")
output_dir <- file.path(user_path, "2022_long_reads", "figure_revision")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

performance <- read.table(
    input_file,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    col.names = c("sample", "depth", "caller", "Sensitivity", "Precision", "F1_score")
)

performance$depth <- as.numeric(performance$depth)
performance$Sensitivity <- as.numeric(performance$Sensitivity)
performance$Precision <- as.numeric(performance$Precision)
performance$F1_score <- as.numeric(performance$F1_score)

caller_levels <- c("TLDR", "xTea", "LOCATE", "TrEMOLO", "GraffiTE")
caller_levels <- c(
    caller_levels[caller_levels %in% performance$caller],
    setdiff(unique(performance$caller), caller_levels)
)
performance$caller <- factor(performance$caller, levels = caller_levels)

metric_data <- rbind(
    data.frame(
        sample = performance$sample,
        depth = performance$depth,
        caller = performance$caller,
        metric = "Sensitivity",
        value = performance$Sensitivity
    ),
    data.frame(
        sample = performance$sample,
        depth = performance$depth,
        caller = performance$caller,
        metric = "Precision",
        value = performance$Precision
    ),
    data.frame(
        sample = performance$sample,
        depth = performance$depth,
        caller = performance$caller,
        metric = "F1-score",
        value = performance$F1_score
    )
)

metric_data$metric <- factor(
    metric_data$metric,
    levels = c("Sensitivity", "Precision", "F1-score")
)

caller_colors <- c(
    TLDR = "#1F77B4",
    xTea = "#2CA02C",
    LOCATE = "#D62728",
    TrEMOLO = "#9467BD",
    GraffiTE = "#FF7F0E"
)
missing_callers <- setdiff(levels(performance$caller), names(caller_colors))
if (length(missing_callers) > 0) {
    extra_colors <- grDevices::hcl.colors(length(missing_callers), palette = "Dark 3")
    names(extra_colors) <- missing_callers
    caller_colors <- c(caller_colors, extra_colors)
}
caller_colors <- caller_colors[levels(performance$caller)]

plot_performance <- function(data) {
    ggplot(
        data,
        aes(x = depth, y = value, color = caller, group = caller)
    ) +
        geom_line(size = 0.45) +
        geom_point(size = 1.6) +
        facet_grid(sample ~ metric) +
        scale_color_manual(values = caller_colors, drop = FALSE) +
        scale_x_continuous(
            breaks = sort(unique(data$depth)),
            expand = expansion(mult = c(0.05, 0.08))
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.2),
            labels = function(x) sprintf("%.1f", x),
            expand = expansion(mult = c(0, 0.03))
        ) +
        labs(
            x = "Sequencing depth",
            y = "Performance",
            color = "Caller"
        ) +
        theme_bw(base_size = 9) +
        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(size = 0.2, color = "grey88"),
            strip.background = element_rect(fill = "grey95", color = "grey75"),
            strip.text = element_text(face = "bold"),
            legend.position = "top",
            legend.title = element_text(face = "bold"),
            legend.key.width = unit(0.7, "cm"),
            axis.text = element_text(color = "black")
        )
}

main_plot <- plot_performance(metric_data)

pdf_file <- file.path(output_dir, paste0(sample, ".performance_by_depth.pdf"))

ggsave(pdf_file, main_plot, width = 7.2, height = 3.2, bg = "white")

message("Saved: ", pdf_file)
