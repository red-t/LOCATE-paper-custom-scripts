
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

result_root <- file.path(user_path, "2022_long_reads", "result")
table_dir <- file.path(user_path, "2022_long_reads", "tables")
figure_dir <- file.path(user_path, "2022_long_reads", "figure_revision")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

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


performance_plot <- function(pre_data, perform, my_legend) {
    pre_data[,"perform"] <- pre_data[,perform]

    caller_list <- as.character(levels(droplevels(pre_data$caller)))

    color_manual <- as.character(color_shape_matrix[caller_list, "color"])
    lines_manual <- as.character(color_shape_matrix[ caller_list, "lines"])
    shape_manual <- as.integer(color_shape_matrix[caller_list, "shape"])
    size_manual <- color_shape_matrix[caller_list, "size"]


    y_max <- max(pre_data[,"perform"])

    if(perform=="F1-score"){
        labels_show <- c("0.0", "0.2", "0.4", "0.6", "0.8","0.9", "" ,"1.0")
        y_lable <- paste0(perform)
    }else{
        labels_show=c(0, 20, 40, 60, 80, 90, "",100)
        
        y_lable <- paste0(perform, " (%)")
    }

    ggplot(pre_data, aes(x = depth, y = perform,  group = caller,  colour = caller, linetype=caller)) +
        geom_point(size = 0.6) +
        geom_line(aes(size=caller)) +
        
        scale_size_manual(values = size_manual) + 
        # scale_shape_manual(values = 1 ) + 
        scale_linetype_manual(values = lines_manual ) + 
        scale_color_manual(values = color_manual)+
        xlab("Sequencing Depth") +
        ylab(y_lable) +
        guides(nrow = 1,size = guide_legend(reverse=TRUE), color = guide_legend(reverse=TRUE), shape = guide_legend(reverse=TRUE), linetype = guide_legend(reverse=TRUE))+
        # guides(nrow = 1, color = guide_legend(reverse=TRUE,label.theme=element_text(size = 1),keywidth=0.5,keyheigh=0.5), shape= guide_legend(reverse=TRUE, keywidth=0.5,keyheigh=0.5,override.aes = list(size = 1),label.theme=element_text(size = 4)))+
        theme(panel.spacing = unit(0.1, "cm")) + 
        theme(legend.position = my_legend, legend.direction = "horizontal") +
        scale_x_continuous(expand = c(0.08,0,0.08,0), breaks = c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), labels = c("1", "", "", "", "", 10, 20, 30, 40, 50)) +
        # scale_y_continuous(breaks = c(0.00, 0.25, 0.5, 0.75, 0.9, 0.95, 1.00), labels = c(0, 25, 50, 75, 90, 95, 100)) +
        scale_y_continuous(expand = c(0.08,0,0.08,0), limits = c(0,1),breaks = c(0.00, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.00), labels = labels_show) +
        # scale_y_continuous(limits = c(0, 1.00),
        #                 breaks = c(0.00, 0.1, 0.2, 0.3, 0.4),
        #                 label = c(0.00, 0.1, 0.2, 0.3, 0.4))+
        # annotate("rect", xmin = 20, xmax = 52, ymin = 0.9, ymax = 1, alpha = 0.3, 
        #      fill = "grey")  +
        mytemp_locate + 
        facet_zoom2(xlim = c(20, 50),ylim=c(0.9,1), split = FALSE,zoom.size = 0.5, show.area=TRUE)
        
        
}

performance_plot_dis <- function(pre_data, perform, my_legend, data_type , TE="All insertion", zoom_in="FALSE") {

    if(TE!="All insertion"){
        pre_data <- pre_data[pre_data$TE==TE, ]
        print(length(pre_data$TE))
        if(length(pre_data$TE) == 0 ){
            return(1)
        }
        plot_title <- TE
    }else{
        plot_title <- ""
    }

    pre_data[,"perform"] <- pre_data[,perform]
    caller_list <- as.character(levels(droplevels(pre_data$caller)))
    
    color_manual <- as.character(color_shape_matrix[caller_list, "color"])
    lines_manual <- as.character(color_shape_matrix[caller_list, "lines"])
    shape_manual <- as.integer(color_shape_matrix[caller_list, "shape"])
    size_manual <- color_shape_matrix[caller_list, "size"]

    y_max <- max(pre_data[,"perform"])

    if(perform == "avg_bp_dis"){
        zoom_ylim <- 5
        y_lable <- paste0("Avg. distance from ", "\n", "detected to breakpoint")
    }

    if(perform == "avg_te_dis"){
        zoom_ylim <- 5
        y_lable <- paste0("Avg. distance from ", "\n", "detected to TE ends")
    }

    if(perform == "avg_div"){
        zoom_ylim <- max(pre_data[pre_data$depth>=20,'avg_div'])
        y_lable <- paste0("Avg. divergency", "\n", "(-log10)")
    }
    p_origin <- ggplot(pre_data, aes(x = depth, y = perform, group = caller , colour = caller, linetype=caller, size=caller)) +
        geom_point(size = 0.6) +
        geom_line(aes(size = caller) ) +
        
        scale_size_manual(values = size_manual) + 
        scale_color_manual(values = color_manual)+
        scale_linetype_manual(values = lines_manual ) + 
        # scale_shape_manual(values = shape_manual) + 
        labs(x="Sequencing Depth", y=y_lable, title = plot_title)+
        guides(nrow = 1, color = guide_legend(reverse=TRUE), shape = guide_legend(reverse=TRUE), linetype = guide_legend(reverse=TRUE))+
        theme(legend.position = my_legend) +
        theme(panel.spacing = unit(0.1, "cm")) + 
        # scale_y_break(breaks = c(5,90),space = 0.2,scales = 1,expand = c(0,0))+
        # geom_hline(aes(yintercept=0), colour="grey", linetype="dashed")+
        scale_y_continuous(expand = c(0.08,0,0.08,0),limits=c(-7,0),breaks = c(-7,-6,-5,-4,-3,-2,-1,0), labels = c(0,-6,-5,-4,-3,-2,-1,0)) +
        # scale_y_continuous(expand = c(0.08,0,0.08,0))+
        scale_x_continuous(expand = c(0.08,0,0.08,0), breaks = c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), labels = c("1", "", "", "", "", 10, 20, 30, 40, 50)) +
        # scale_y_continuous(limits = c(0, y_max), breaks = c(0.00,0.001, 0.002, 0.005, 0.01, 0.025), labels = c(0.00,0.001, 0.002, 0.005, 0.01, 0.025)) +
        mytemp_locate
    # p_origin <- gg.gap(plot = p_origin.1,
    #        segments = c(5, 99),
    #        tick_width = c(1,1),
    #        ylim = c(0, 100)) + mytemp_locate

    if(zoom_in=="TRUE"){
        p_return <- p_origin + facet_zoom2(xlim = c(20, 50), ylim=c(0,zoom_ylim), split = FALSE, zoom.size = 0.5)
    }else{
        p_return <- p_origin
    }
    return(p_return)
}

### perfomance

genome <- "hg38"
data_path <- file.path(result_root, "simulation", genome)


perfomance <- read.table(paste0(data_path,"/intersection_revision/pp_intersection/data/performance.txt"))
# perfomance <- read.table(paste0(data_path,"/intersection/pp_intersection/performance_ex_internal.txt"))
colnames(perfomance) <- c("data_type", "depth", "caller", "Sensitivity", "Precision", "F1-score")

perfomance_d.sensitivity <- dcast(perfomance, data_type + depth~caller, value.var = "Sensitivity" ) # c("Sensitivity", "Precision", "F1-score")
perfomance_d.sensitivity$performance <- "Sensitivity"
perfomance_d.precision <- dcast(perfomance, data_type + depth~caller, value.var = "Precision" ) 
perfomance_d.precision$performance <- "Precision"
perfomance_d.f1_score <- dcast(perfomance, data_type + depth~caller, value.var = "F1-score" )
perfomance_d.f1_score$performance <- "F1-score"

perfomance_d <- rbind(perfomance_d.sensitivity, perfomance_d.precision, perfomance_d.f1_score )  

write.xlsx(perfomance_d, file.path(table_dir, "hg38.performance_metrics.xlsx"), sheetName="hg38.performance", colNames = TRUE)

perfomance <- perfomance[ ! ( perfomance$caller=="TELR" & perfomance$Precision=='0' ), ]
perfomance <- perfomance[ ! ( perfomance$caller=="MELT" & perfomance$Precision=='0' ), ]
perfomance <- perfomance[ ! ( perfomance$caller=="TrEMOLO" & perfomance$Precision=='0' ), ]
data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))


# caller_hg38 <- c("TEMP2", "MELT","TLDR_all","TLDR_pass" , "xTea_all","xTea_pass",  "TELR", "TrEMOLO",  "MEHunter",  "GraffiTE", "LOCATE_all", "LOCATE_pass")

caller_hg38 <- c("TEMP2", "MELT","TLDR_pass" , "xTea_pass",  "TELR", "TrEMOLO", "PALMER",  "MEHunter",   "GraffiTE", "LOCATE_pass")
perfomance <- perfomance[ perfomance$caller!="TLDR_all" & perfomance$caller!="xTea_all" & perfomance$caller!="LOCATE_all", ]
perfomance$caller <- factor(perfomance$caller, levels=caller_hg38)



performance_plot_spe <- function(data_type){
    if (data_type == "simulation_germ_ccs"){
        fig_title <- "Simulated Pacbio CCS"
    }
    if (data_type == "simulation_germ_clr"){
        fig_title <- "Simulated Pacbio CLR"
    }
    if (data_type == "simulation_germ_ont"){
        fig_title <- "Simulated Nanopore"
    }
    pre_data <- perfomance[perfomance$data_type==data_type,]                      
    p1 <- performance_plot(pre_data, "Sensitivity","None")
    p2 <- performance_plot(pre_data, "Precision","None")
    p3 <- performance_plot(pre_data, "F1-score","None")
    p_legend <- performance_plot(pre_data, "Sensitivity","top")
    legend <- get_plot_component(p_legend, 'guide-box-top', return_all = TRUE)

    pic <- cowplot::plot_grid(p1, p2, p3,nrow = 1, ncol=3,labels = c("A", "B", "C"))
    title <- ggdraw() + draw_label(paste0(fig_title), size=10)
    plot_grid(title, pic, legend, ncol=1, rel_heights=c(0.3, 1.2, 0.6),scale = c(1, 1, 1)) # rel_heights values control title margins

}

data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))
pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_spe(x[1]) )
pic <- cowplot::plot_grid(pic_all[[1]], pic_all[[2]], pic_all[[3]],nrow = 3, ncol=1, label_size = 18)

ggsave(pic, filename = file.path(figure_dir, "hg38.performance_by_depth.pdf"), width = 7.5, height = 8, bg = 'transparent')




# read data
genome <- "hg38"
data_path <- file.path(result_root, "simulation", genome)
# data_path <- file.path(result_root, "simulation", genome)
perfomance_dis <- read.table(paste0(data_path, "/intersection_revision/all.merge.dis.txt"), header = T)
perfomance_dis <- perfomance_dis[perfomance_dis$frequency>0.2,]
head(perfomance_dis)


### calculated
# bp_dis
perfomance.bp_dis <- melt(perfomance_dis, id.vars=c("insID", "sample", "depth"), measure.vars =c("LOCATE_all.bp_dis", "LOCATE_pass.bp_dis", "MELT.bp_dis","TEMP2.bp_dis", "TLDR_all.bp_dis","TLDR_pass.bp_dis" , "xTea_all.bp_dis", "xTea_pass.bp_dis", "PALMER.bp_dis", "TELR.bp_dis", "TrEMOLO.bp_dis", "MEHunter.bp_dis", "GraffiTE.bp_dis"), variable.name = "caller", value.name = "bp_dis")
perfomance.bp_dis$caller <- sub("\\..*", "", perfomance.bp_dis$caller)
perfomance.bp_dis <- perfomance.bp_dis[perfomance.bp_dis$bp_dis != "-", ]
perfomance.bp_dis$bp_dis <- as.numeric(perfomance.bp_dis$bp_dis)

# te_dis
perfomance.te_dis <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.te_dis", "LOCATE_pass.te_dis", "MELT.te_dis","TEMP2.te_dis", "TLDR_all.te_dis","TLDR_pass.te_dis" , "xTea_all.te_dis", "xTea_pass.te_dis", "PALMER.te_dis", "TELR.te_dis", "TrEMOLO.te_dis", "MEHunter.te_dis", "GraffiTE.te_dis" ), variable.name = "caller", value.name = "te_dis")
perfomance.te_dis$caller <- sub("\\..*", "", perfomance.te_dis$caller)
perfomance.te_dis <- perfomance.te_dis[perfomance.te_dis$te_dis != "-", ]
perfomance.te_dis$te_dis <- as.numeric(perfomance.te_dis$te_dis)

# te_identity
perfomance.te_identity <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.te_identity", "LOCATE_pass.te_identity", "MELT.te_identity","TEMP2.te_identity", "TLDR_all.te_identity","TLDR_pass.te_identity" , "xTea_all.te_identity", "xTea_pass.te_identity", "PALMER.te_identity", "TELR.te_identity", "TrEMOLO.te_identity", "MEHunter.te_identity", "GraffiTE.te_identity" ), variable.name = "caller", value.name = "te_identity")
perfomance.te_identity$caller <- sub("\\..*", "", perfomance.te_identity$caller)
perfomance.te_identity <- perfomance.te_identity[perfomance.te_identity$te_identity != "-", ]
perfomance.te_identity$te_identity <- as.numeric(perfomance.te_identity$te_identity)

# genotype
perfomance.genotype <- melt(perfomance_dis, id.vars=c("insID", "sample", "depth"), measure.vars =c("LOCATE_all.genotype", "LOCATE_pass.genotype", "MELT.genotype","TEMP2.genotype", "TLDR_all.genotype","TLDR_pass.genotype" , "xTea_all.genotype", "xTea_pass.genotype", "PALMER.genotype", "TELR.genotype", "TrEMOLO.genotype", "MEHunter.genotype", "GraffiTE.genotype"), variable.name = "caller", value.name = "genotype")
perfomance.genotype$caller <- sub("\\..*", "", perfomance.genotype$caller)
perfomance.genotype <- perfomance.genotype[perfomance.genotype$genotype != "-", ]
perfomance.genotype$genotype <- as.numeric(perfomance.genotype$genotype)


# div
perfomance.te_len <- melt(perfomance_dis, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_all.te_len", "LOCATE_pass.te_len", "MELT.te_len","TEMP2.te_len", "TLDR_all.te_len","TLDR_pass.te_len" , "xTea_all.te_len", "xTea_pass.te_len", "PALMER.te_len", "TELR.te_len", "TrEMOLO.te_len", "MEHunter.te_len", "GraffiTE.te_len" ), variable.name = "caller", value.name = "te_len")
perfomance.te_len$caller <- sub("\\..*", "", perfomance.te_len$caller)
perfomance.div <- melt(perfomance_dis, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_all.div", "LOCATE_pass.div", "MELT.div","TEMP2.div", "TLDR_all.div","TLDR_pass.div" , "xTea_all.div", "xTea_pass.div", "PALMER.div", "TELR.div", "TrEMOLO.div", "MEHunter.div" , "GraffiTE.div"), variable.name = "caller", value.name = "div")
perfomance.div$caller <- sub("\\..*", "", perfomance.div$caller)
perfomance.len_div <- merge(perfomance.te_len, perfomance.div, by = c("insID","sample", "depth", "caller", "TE"))
perfomance.len_div <- perfomance.len_div[perfomance.len_div$div != "-", ]
perfomance.len_div$div <- as.numeric(perfomance.len_div$div)
perfomance.len_div$te_len <- as.numeric(perfomance.len_div$te_len)

# LOCATE_TLDR
perfomance_dis.LOCATE_TLDR <- perfomance_dis[perfomance_dis$LOCATE_pass.div!="-" & perfomance_dis$TLDR_pass.div!="-", ]
perfomance.LOCATE_TLDR.te_len <- melt(perfomance_dis.LOCATE_TLDR, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_pass.te_len","TLDR_pass.te_len"  ), variable.name = "caller", value.name = "te_len")
perfomance.LOCATE_TLDR.te_len$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR.te_len$caller)
perfomance.LOCATE_TLDR.div <- melt(perfomance_dis.LOCATE_TLDR, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_pass.div","TLDR_pass.div" ), variable.name = "caller", value.name = "div")
perfomance.LOCATE_TLDR.div$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR.div$caller)
perfomance.LOCATE_TLDR.len_div <- merge(perfomance.LOCATE_TLDR.te_len, perfomance.LOCATE_TLDR.div, by = c("insID","sample", "depth", "caller", "TE"))
perfomance.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[perfomance.LOCATE_TLDR.len_div$div != "-", ]
perfomance.LOCATE_TLDR.len_div$div <- as.numeric(perfomance.LOCATE_TLDR.len_div$div)
perfomance.LOCATE_TLDR.len_div$te_len <- as.numeric(perfomance.LOCATE_TLDR.len_div$te_len)

# LOCATE_TLDR_xTea
perfomance_dis.LOCATE_TLDR_xTea <- perfomance_dis[perfomance_dis$LOCATE_pass.div!="-" & perfomance_dis$TLDR_pass.div!="-" & perfomance_dis$xTea_pass.div!="-", , ]
perfomance.LOCATE_TLDR_xTea.te_len <- melt(perfomance_dis.LOCATE_TLDR_xTea, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_pass.te_len","TLDR_pass.te_len" ,"xTea_pass.te_len"  ), variable.name = "caller", value.name = "te_len")
perfomance.LOCATE_TLDR_xTea.te_len$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_xTea.te_len$caller)
perfomance.LOCATE_TLDR_xTea.div <- melt(perfomance_dis.LOCATE_TLDR_xTea, id.vars=c("insID","sample", "depth", "TE"), measure.vars =c("LOCATE_pass.div","TLDR_pass.div","xTea_pass.div" ), variable.name = "caller", value.name = "div")
perfomance.LOCATE_TLDR_xTea.div$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_xTea.div$caller)
perfomance.LOCATE_TLDR_xTea.len_div <- merge(perfomance.LOCATE_TLDR_xTea.te_len, perfomance.LOCATE_TLDR_xTea.div, by = c("insID","sample", "depth", "caller", "TE"))
perfomance.LOCATE_TLDR_xTea.len_div <- perfomance.LOCATE_TLDR_xTea.len_div[perfomance.LOCATE_TLDR_xTea.len_div$div != "-", ]
perfomance.LOCATE_TLDR_xTea.len_div$div <- as.numeric(perfomance.LOCATE_TLDR_xTea.len_div$div)
perfomance.LOCATE_TLDR_xTea.len_div$te_len <- as.numeric(perfomance.LOCATE_TLDR_xTea.len_div$te_len)

perfomance.bp_dis <- perfomance.bp_dis[! perfomance.bp_dis$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.te_dis <- perfomance.te_dis[! perfomance.te_dis$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.te_identity <- perfomance.te_identity[! perfomance.te_identity$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.len_div <- perfomance.len_div[! perfomance.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[! perfomance.LOCATE_TLDR.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.LOCATE_TLDR_xTea.len_div <- perfomance.LOCATE_TLDR_xTea.len_div[! perfomance.LOCATE_TLDR_xTea.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]


# levels
caller_hg38_dis <- c("TEMP2", "MELT","TLDR_pass" , "xTea_pass",  "TELR", "TrEMOLO", "PALMER", "MEHunter", "GraffiTE","LOCATE_pass")
perfomance.bp_dis$caller <- factor(perfomance.bp_dis$caller, levels=caller_hg38_dis)
perfomance.te_dis$caller <- factor(perfomance.te_dis$caller, levels=caller_hg38_dis)
perfomance.te_identity$caller <- factor(perfomance.te_identity$caller, levels=caller_hg38_dis)
perfomance.len_div$caller <- factor(perfomance.len_div$caller, levels=caller_hg38_dis)
perfomance.LOCATE_TLDR.len_div$caller <- factor(perfomance.LOCATE_TLDR.len_div$caller, levels=caller_hg38_dis)
perfomance.LOCATE_TLDR_xTea.len_div$caller <- factor(perfomance.LOCATE_TLDR_xTea.len_div$caller, levels=caller_hg38_dis)

# perfomance.LOCATE_TLDR_xTea.len_div.1 <- dcast(perfomance.LOCATE_TLDR_xTea.len_div, data_type + depth~caller, value.var = "F1-score" )


# plot divergency

performance_plot_dis_spe <- function(data_type){
    pre_data.bp_dis <- perfomance.bp_dis[perfomance.bp_dis$sample==data_type,]
    pre_data.bp_dis <- pre_data.bp_dis %>% group_by(caller, depth) %>% dplyr::summarise(avg_bp_dis = mean(bp_dis), .groups = "keep" )

    pre_data.te_identity <- perfomance.te_identity[perfomance.te_identity$sample==data_type,]
    pre_data.te_identity <- pre_data.te_identity %>% group_by(caller, depth) %>% dplyr::summarise(avg_te_identity = mean(te_identity), .groups = "keep")

    pre_data.len_div <- perfomance.len_div[perfomance.len_div$sample==data_type,]

    pre_data.len_div <- pre_data.len_div %>% group_by(caller, depth) %>% dplyr::summarise(avg_div = sum(te_len*div)/sum(te_len), .groups = "keep")
    
    pre_data.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[perfomance.LOCATE_TLDR.len_div$sample==data_type,]
    pre_data.LOCATE_TLDR.len_div <- pre_data.LOCATE_TLDR.len_div %>% group_by(caller, depth) %>% dplyr::summarise(avg_div = sum(te_len*div)/sum(te_len), .groups = "keep")
    

    pre_data.LOCATE_TLDR_xTea.len_div <- perfomance.LOCATE_TLDR_xTea.len_div[perfomance.LOCATE_TLDR_xTea.len_div$sample==data_type,]
    pre_data.LOCATE_TLDR_xTea.len_div <- pre_data.LOCATE_TLDR_xTea.len_div %>% group_by(caller, depth) %>% dplyr::summarise(avg_div = sum(te_len*div)/sum(te_len), .groups = "keep")
    pre_data.LOCATE_TLDR_xTea.len_div.1 <-  dcast( pre_data.LOCATE_TLDR_xTea.len_div, depth~caller, value.var = "avg_div" )
    # pre_data.LOCATE_TLDR_xTea.len_div$avg_div <- -log(pre_data.LOCATE_TLDR_xTea.len_div$avg_div + 0.000001, 10)
    pre_data.LOCATE_TLDR_xTea.len_div <- pre_data.LOCATE_TLDR_xTea.len_div %>% mutate(avg_div=case_when(avg_div==0 ~ -7,
                                                                                                        TRUE ~ log(avg_div,10)))  # <- -log(pre_data.LOCATE_TLDR_xTea.len_div$avg_div + 0.000001, 10)
    
    # write.table(pre_data.LOCATE_TLDR_xTea.len_div, paste0(data_path, "/figure/",genome,".",data_type, ".div.txt"),row.names = FALSE, sep ="\t",col.names =FALSE, quote =FALSE)
    # pre_data.LOCATE_TLDR_xTea.len_div$avg_div[is.na(pre_data.LOCATE_TLDR_xTea.len_div$avg_div)] <- "NA"
    
    
    # sheet <- createSheet(wb, sheetName = paste0(data_type,".hg38.div"))
    # writeData(wb, sheet, pre_data.LOCATE_TLDR_xTea.len_div)


    # p3 <- performance_plot_dis(pre_data.len_div, "avg_div","None",data_type,"All insertion" , "FALSE")
    # p4 <- performance_plot_dis(pre_data.LOCATE_TLDR.len_div, "avg_div","None",data_type)
    print(table((pre_data.LOCATE_TLDR_xTea.len_div$caller)))
    p5 <- performance_plot_dis(pre_data.LOCATE_TLDR_xTea.len_div, "avg_div","None",data_type,"All insertion" ,"FALSE")

    p_legend <- performance_plot_dis(pre_data.LOCATE_TLDR_xTea.len_div, "avg_div","top",data_type)
    legend <- get_plot_component(p_legend, 'guide-box-top', return_all = TRUE)
    

    # pic <- cowplot::plot_grid(p3, p5, nrow = 1, ncol=2, rel_widths = c(1.15,1))
    # pic <- cowplot::plot_grid(p5, nrow = 1, ncol=1)

    if (data_type == "simulation_germ_ccs"){
        fig_title <- "Simulated Pacbio HiFi"
    }
    if (data_type == "simulation_germ_clr"){
        fig_title <- "Simulated Pacbio CLR"
    }
    if (data_type == "simulation_germ_ont"){
        fig_title <- "Simulated Nanopore"
    }

    title <- ggdraw() + draw_label(paste0(fig_title), size=10)
    p.temp <- plot_grid(title, p5, ncol=1, rel_heights=c(0.3, 1.2))
    # return(p.temp)
    return(list(p.temp,legend ,pre_data.LOCATE_TLDR_xTea.len_div.1))# rel_heights values control title margins

}


# wb <- createWorkbook()

data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))
pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_dis_spe(x[1]) )
pic <- cowplot::plot_grid(pic_all[[1]][[1]], pic_all[[2]][[1]],pic_all[[3]][[1]],nrow = 1, ncol=3)
# pic.1 <- cowplot::plot_grid(pic, pic_all[[1]][[2]],nrow = 2, ncol=1, rel_heights=c(1, 0.4),scale = c(1,1))

ggsave(pic, filename = file.path(figure_dir, "hg38.sequence_divergence_by_depth.pdf"), width = 6, height = 2.4)

## genotype

perfomance.genotype <- melt(perfomance_dis, id.vars=c("insID", "sample", "depth"), measure.vars =c("LOCATE_all.genotype", "LOCATE_pass.genotype", "MELT.genotype","TEMP2.genotype", "TLDR_all.genotype","TLDR_pass.genotype" , "xTea_all.genotype", "xTea_pass.genotype", "PALMER.genotype", "TELR.genotype", "TrEMOLO.genotype", "MEHunter.genotype", "GraffiTE.genotype"), variable.name = "caller", value.name = "genotype")
perfomance.genotype$caller <- sub("\\..*", "", perfomance.genotype$caller)
perfomance.genotype <- perfomance.genotype[perfomance.genotype$genotype != "-", ]
perfomance.genotype$genotype <- as.numeric(perfomance.genotype$genotype)


                 #!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

output_dir <- figure_dir
depth_levels <- c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50)

# script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
# if (is.na(script_path)) {
#   args <- commandArgs(trailingOnly = FALSE)
#   file_arg <- grep("^--file=", args, value = TRUE)
#   script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else NA_character_
# }

# genotype_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
genotype_dir <- file.path(result_root, "simulation", "hg38", "intersection_revision")
count_file_candidates <- unique(c(
  file.path(genotype_dir, "all.merge.count.txt"),
  file.path(getwd(), "all.merge.count.txt"),
  file.path(getwd(), "genotype", "all.merge.count.txt")
))
                        
default_count_file <- count_file_candidates[file.exists(count_file_candidates)][1]
if (is.na(default_count_file)) {
  default_count_file <- count_file_candidates[1]
}

callers <- c(
  "TEMP2", "MELT", "TLDR_pass", "TLDR_all", "xTea_pass", "xTea_all",
  "TELR", "TrEMOLO", "PALMER", "MEHunter", "GraffiTE",
  "LOCATE_all", "LOCATE_pass"
)

colors <- c(
  "#FFC107", "#FF6F00", "#2196F3", "#2196F3", "#4CAF50", "#4CAF50",
  "#BDBDBD", "#8D6E63", "#3F51B5", "#7E57C2", "#f1b6da",
  "red", "red"
)

line_types <- c(
  "solid", "solid", "solid", "dashed", "solid", "dashed",
  "solid", "solid", "solid", "solid", "solid",
  "dashed", "solid"
)

line_size <- c(
  0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
  0.2, 0.2, 0.2, 0.2, 0.2,
  0.4, 0.4
)

color_shape_matrix <- data.frame(
  caller = callers,
  color = colors,
  lines = line_types,
  size = line_size,
  stringsAsFactors = FALSE
)
rownames(color_shape_matrix) <- callers

genotype_accuracy_data <- NULL

normalize_caller_name <- function(caller) {
  gsub("-", "_", as.character(caller))
}

sample_label <- function(sample_name) {
  sample_name <- as.character(sample_name)
  labels <- c(
    simulation_germ_ccs = "HiFi",
    simulation_germ_clr = "CLR",
    simulation_germ_ont = "ONT"
  )
  label <- labels[sample_name]
  unname(ifelse(is.na(label), sample_name, label))
}

format_count_keys <- function(rows, max_rows = 10) {
  rows <- rows[seq_len(min(nrow(rows), max_rows)), , drop = FALSE]
  keys <- paste0(rows$sample, "/depth=", rows$depth, "/caller=", rows$caller)
  paste(keys, collapse = "; ")
}

read_count_data <- function(count_file = default_count_file) {
  if (!file.exists(count_file)) {
    stop("Count file not found: ", count_file)
  }

  read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
}

prepare_count_data <- function(count_df) {
  if (!"meth" %in% names(count_df) && "caller" %in% names(count_df)) {
    count_df$meth <- count_df$caller
  }

  required_cols <- c("sample", "depth", "meth", "caller_count", "gold_count")
  missing_cols <- setdiff(required_cols, names(count_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in count data: ", paste(missing_cols, collapse = ", "))
  }

  count_data <- count_df %>%
    dplyr::mutate(
      sample = as.character(sample),
      caller = normalize_caller_name(meth),
      depth = as.integer(depth),
      caller_count = as.integer(caller_count),
      gold_count = as.integer(gold_count)
    ) %>%
    dplyr::filter(depth %in% depth_levels) %>%
    dplyr::select(sample, caller, depth, caller_count, gold_count)

  duplicate_rows <- count_data %>%
    dplyr::count(sample, caller, depth, name = "rows") %>%
    dplyr::filter(rows > 1)

  if (nrow(duplicate_rows) > 0) {
    stop("Duplicate count rows found: ", format_count_keys(duplicate_rows))
  }

  count_data
}

add_union_accuracy <- function(summary_data, count_data) {
  accuracy_data <- summary_data %>%
    dplyr::left_join(count_data, by = c("sample", "caller", "depth"))

  missing_rows <- accuracy_data %>%
    dplyr::filter(is.na(caller_count) | is.na(gold_count))

  if (nrow(missing_rows) > 0) {
    stop("Missing caller/gold counts for: ", format_count_keys(missing_rows))
  }

  invalid_rows <- accuracy_data %>%
    dplyr::filter(correct > caller_count | correct > gold_count)

  if (nrow(invalid_rows) > 0) {
    warning(
      "Correct count exceeds caller_count or gold_count for: ",
      format_count_keys(invalid_rows),
      call. = FALSE
    )
  }

  accuracy_data %>%
    dplyr::mutate(
      union_count = caller_count + gold_count - n,
      accuracy = dplyr::if_else(union_count > 0, correct / union_count, NA_real_),
      no_call_rate = dplyr::if_else(n > 0, no_call / n, NA_real_)
    ) %>%
    dplyr::select(
      sample, caller, depth, n, correct, wrong, no_call,
      caller_count, gold_count, union_count, accuracy, no_call_rate
    )
}

prepare_genotype_accuracy <- function(genotype_df, count_df = NULL, count_file = default_count_file) {
  required_cols <- c("sample", "depth", "caller", "genotype")
  missing_cols <- setdiff(required_cols, names(genotype_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(count_df)) {
    count_df <- read_count_data(count_file)
  }
  count_data <- prepare_count_data(count_df)

  genotype_df %>%
    dplyr::mutate(
      sample = as.character(sample),
      caller = normalize_caller_name(caller),
      depth = as.integer(depth),
      genotype = as.integer(genotype)
    ) %>%
    dplyr::filter(depth %in% depth_levels) %>%
    dplyr::group_by(sample, caller, depth) %>%
    dplyr::summarise(
      n = dplyr::n(),
      correct = sum(genotype == 1, na.rm = TRUE),
      wrong = sum(genotype == 0, na.rm = TRUE),
      no_call = sum(genotype == 2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    add_union_accuracy(count_data) %>%
    dplyr::arrange(sample, caller, depth)
}

set_genotype_accuracy_data <- function(genotype_df, count_df = NULL, count_file = default_count_file) {
  genotype_accuracy_data <<- prepare_genotype_accuracy(
    genotype_df,
    count_df = count_df,
    count_file = count_file
  )
  invisible(genotype_accuracy_data)
}

performance_plot_spe <- function(sample_name, my_legend = "bottom") {
  if (is.null(genotype_accuracy_data)) {
    stop("Run set_genotype_accuracy_data(genotype_df) before calling performance_plot_spe().")
  }
  if (!exists("mytemp_locate", inherits = TRUE)) {
    stop("mytemp_locate is not found. Define your theme object before calling performance_plot_spe().")
  }

  plot_data <- genotype_accuracy_data %>%
    dplyr::filter(sample == sample_name)

  if (nrow(plot_data) == 0) {
    stop("No data found for sample: ", sample_name)
  }

  caller_list <- intersect(rownames(color_shape_matrix), unique(as.character(plot_data$caller)))

  plot_data <- plot_data %>%
    dplyr::filter(caller %in% caller_list) %>%
    dplyr::mutate(caller = factor(caller, levels = caller_list))

  color_manual <- as.character(color_shape_matrix[caller_list, "color"])
  lines_manual <- as.character(color_shape_matrix[caller_list, "lines"])
  size_manual <- color_shape_matrix[caller_list, "size"]

  p <- ggplot(
    plot_data,
    aes(x = depth, y = accuracy, group = caller, colour = caller, linetype = caller)
  ) +
    geom_point(size = 0.6) +
    geom_line(aes(size = caller)) +
    scale_size_manual(values = size_manual) +
    scale_linetype_manual(values = lines_manual) +
    scale_color_manual(values = color_manual) +
    scale_x_continuous(
      expand = c(0.08, 0, 0.08, 0),
      breaks = depth_levels,
      labels = c("1", "", "", "", "", "10", "20", "30", "40", "50")
    ) +
    scale_y_continuous(
      expand = c(0.08, 0, 0.08, 0),
      limits = c(0, 1),
      breaks = c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0),
      labels = c(0, 20, 40, 60, 80, 90, "", 100)
    ) +
    labs(
      title = sample_label(sample_name),
      x = "Sequencing Depth",
      y = "Genotype accuracy (%)"
    ) +
    guides(
      size = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE),
      colour = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE),
      linetype = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)
    ) +
    theme(panel.spacing = unit(0.1, "cm")) +
    theme(legend.position = my_legend, legend.direction = "horizontal") +
    mytemp_locate

  if (requireNamespace("ggforce", quietly = TRUE)) {
    p <- p + ggforce::facet_zoom(
      xlim = c(20, 50),
      ylim = c(0.9, 1),
      split = FALSE,
      zoom.size = 0.5,
      show.area = TRUE
    )
  }

  p
}

plot_genotype_accuracy <- function(
  genotype_df,
  output_prefix = "hg38.genotype_accuracy",
  count_df = NULL,
  count_file = default_count_file
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  accuracy_data <- set_genotype_accuracy_data(
    genotype_df,
    count_df = count_df,
    count_file = count_file
  )

  write.csv(
    accuracy_data,
    file.path(output_dir, paste0(output_prefix, ".summary.csv")),
    row.names = FALSE
  )

  for (sample_name in unique(accuracy_data$sample)) {
    p <- performance_plot_spe(sample_name)

    ggsave(
      filename = file.path(output_dir, paste0(output_prefix, ".", sample_name, ".pdf")),
      plot = p,
      width = 4.2,
      height = 3.0,
      units = "in",
      device = cairo_pdf
    )

    ggsave(
      filename = file.path(output_dir, paste0(output_prefix, ".", sample_name, ".png")),
      plot = p,
      width = 4.2,
      height = 3.0,
      units = "in",
      dpi = 300
    )
  }

  if (requireNamespace("cowplot", quietly = TRUE)) {
    sample_order <- c("simulation_germ_ccs", "simulation_germ_clr", "simulation_germ_ont")
    sample_order <- intersect(sample_order, unique(accuracy_data$sample))

    if (length(sample_order) > 0) {
      shared_legend <- cowplot::get_legend(performance_plot_spe(sample_order[1], my_legend = "bottom"))
      pic_all <- lapply(sample_order, function(sample_name) {
        performance_plot_spe(sample_name, my_legend = "none")
      })
      pic_row <- cowplot::plot_grid(plotlist = pic_all, nrow = 1, ncol = length(pic_all))
      pic <- cowplot::plot_grid(pic_row, shared_legend, ncol = 1, rel_heights = c(1, 0.14))

      ggsave(
        filename = file.path(output_dir, paste0(output_prefix, ".combined.pdf")),
        plot = pic,
        width = 12.6,
        height = 3.4,
        units = "in",
        device = cairo_pdf
      )

      ggsave(
        filename = file.path(output_dir, paste0(output_prefix, ".combined.png")),
        plot = pic,
        width = 12.6,
        height = 3.4,
        units = "in",
        dpi = 300
      )
    }
  }

  invisible(accuracy_data)
}

# Usage:
# genotype_df must contain at least these columns:
# sample, depth, caller, genotype
# count_df/count_file must contain:
# sample, depth, meth (or caller), caller_count, gold_count
#
# Define mytemp_locate before plotting.
#
# Example:
# genotype_df <- read.csv("your_genotype_result.csv")
# set_genotype_accuracy_data(genotype_df, count_file = "all.merge.count.txt")
# sample_order <- c("simulation_germ_ccs", "simulation_germ_clr", "simulation_germ_ont")
# shared_legend <- cowplot::get_legend(performance_plot_spe(sample_order[1], my_legend = "bottom"))
# pic_all <- lapply(sample_order, function(x) performance_plot_spe(x, my_legend = "none"))
# pic <- cowplot::plot_grid(cowplot::plot_grid(plotlist = pic_all, nrow = 1), shared_legend, ncol = 1)
#
# Or save summary, separate plots, and combined plot:
# plot_genotype_accuracy(genotype_df, count_file = "all.merge.count.txt")

   
plot_genotype_accuracy(perfomance.genotype)
set_genotype_accuracy_data(perfomance.genotype)


data_type_matrix <- data.frame(data_type = c("simulation_germ_ccs", "simulation_germ_clr", "simulation_germ_ont"))
pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_spe(x[1]))
legend <- get_plot_component(pic_all[[1]], 'guide-box-bottom', return_all = TRUE)

pic <- cowplot::plot_grid(pic_all[[1]], pic_all[[2]], pic_all[[3]], nrow = 1, ncol = 3, label_size = 18)

                 
ggsave(pic, filename = file.path(figure_dir, "hg38.genotype_accuracy.pdf"), width = 9, height = 3, bg = 'transparent')
              

                 
                 
