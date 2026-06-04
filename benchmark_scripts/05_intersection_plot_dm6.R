
# detach("package:magrittr", unload = TRUE)
library("ggplot2")
library("cowplot")
library("dplyr")
library("reshape2")
library("ggbreak")
library("xlsx")
library("openxlsx")

user_path <- Sys.getenv("user_path", unset = Sys.getenv("USER_PATH", unset = "/path/to/user"))
# install.packages('magrittr')
# library("magrittr");
# install.packages("dplyr")
library("dplyr")
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


performance_plot <- function(pre_data, perform, my_legend) {
    pre_data[,"perform"] <- pre_data[,perform]

    caller_list <- as.character(levels(droplevels(pre_data$caller)))

    color_manual <- as.character(color_shape_matrix[caller_list, "color"])
    lines_manual <- as.character(color_shape_matrix[ caller_list, "lines"])
    shape_manual <- as.integer(color_shape_matrix[caller_list, "shape"])
    size_manual <- color_shape_matrix[caller_list, "size"]

    print(size_manual)

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
        facet_zoom2(xlim = c(20, 50),ylim=c(0.8,1), split = FALSE,zoom.size = 0.5, show.area=TRUE)
        
        
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
        # scale_y_continuous(expand = c(0.08,0,0.08,0), limits=c(1,8),breaks = c(1,2,3,4,5,6,7), labels = c(1,2,3,4,5,6,"X")) +
        # scale_y_continuous(expand = c(0.08,0,0.08,0),limits=c(-7,0),breaks = c(-7,-6,-5,-4,-3,-2,-1,0), labels = c(0,-6,-5,-4,-3,-2,-1,0)) +

        scale_y_continuous(expand = c(0.08,0,0.08,0), limits=c(0,30),breaks = c(0,1,2,3,4,5,10,20,30), labels = c(0,1,2,3,4,5,10,20,30))+
        scale_x_continuous(expand = c(0.08,0,0.08,0), breaks = c(1, 2, 3, 4, 5, 10, 20, 30, 40, 50), labels = c("1", "", "", "", "", 10, 20, 30, 40, 50)) +
        # scale_y_continuous(limits = c(0, y_max), breaks = c(0.00,0.001, 0.002, 0.005, 0.01, 0.025), labels = c(0.00,0.001, 0.002, 0.005, 0.01, 0.025)) +
        mytemp_locate
    # p_origin <- gg.gap(plot = p_origin.1,
    #        segments = c(5, 99),
    #        tick_width = c(1,1),
    #        rel_heights = c(1, 0.001, 0.2),    #        ylim = c(0, 100)) + mytemp_locate

    if(zoom_in=="TRUE"){
        p_return <- p_origin + facet_zoom2(xlim = c(20, 50), ylim=c(0,zoom_ylim), split = FALSE, zoom.size = 0.5)
    }else{
        p_return <- p_origin
    }
    return(p_return)
}



### plot performance

genome <- "dm6"
result_root <- file.path(user_path, "2022_long_reads", "result")
table_dir <- file.path(user_path, "2022_long_reads", "tables")
figure_dir <- file.path(user_path, "2022_long_reads", "figure")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
data_path <- file.path(result_root, "simulation", genome)

perfomance <- read.table(paste0(data_path,"/intersection/pp_intersection/data/performance.txt"))
colnames(perfomance) <- c("data_type", "depth", "caller", "Sensitivity", "Precision", "F1-score")
perfomance <- perfomance[ perfomance$caller!="TLDR_all" & perfomance$caller!="LOCATE_all", ]

perfomance_d.sensitivity <- dcast(perfomance, data_type + depth~caller, value.var = "Sensitivity" ) # c("Sensitivity", "Precision", "F1-score")
perfomance_d.sensitivity$performance <- "Sensitivity"
perfomance_d.precision <- dcast(perfomance, data_type + depth~caller, value.var = "Precision" ) 
perfomance_d.precision$performance <- "Precision"
perfomance_d.f1_score <- dcast(perfomance, data_type + depth~caller, value.var = "F1-score" )
perfomance_d.f1_score$performance <- "F1-score"

perfomance_d <- rbind(perfomance_d.sensitivity, perfomance_d.precision, perfomance_d.f1_score )  

write.xlsx(perfomance_d, file.path(table_dir, "dm6.performance_metrics.xlsx"), sheetName="dm6.performance", colNames = TRUE)

caller_dm6_dis <- c("TEMP2", "MELT","TLDR_pass","TELR", "TrEMOLO", "GraffiTE","LOCATE_pass")

perfomance$caller <- factor(perfomance$caller, levels=caller_dm6_dis)
perfomance <- perfomance[ ! ( perfomance$caller=="TELR" & perfomance$Precision=='0' ), ]
perfomance <- perfomance[ ! ( perfomance$caller=="MELT" & perfomance$Precision=='0' ), ]

performance_plot_spe <- function(data_type){
    pre_data <- perfomance[perfomance$data_type==data_type,]                      
    p1 <- performance_plot(pre_data, "Sensitivity","None")
    p2 <- performance_plot(pre_data, "Precision","None")
    p3 <- performance_plot(pre_data, "F1-score","None")
    p_legend <- performance_plot(pre_data, "Sensitivity","top")
    legend <- get_plot_component(p_legend, 'guide-box-top', return_all = TRUE)
    if (data_type == "simulation_germ_ccs"){
        fig_title <- "HiFi"
    }
    if (data_type == "simulation_germ_clr"){
        fig_title <- "CLR"
    }
    if (data_type == "simulation_germ_ont"){
        fig_title <- "ONT"
    }
    pic <- cowplot::plot_grid(p1, p2, p3,nrow = 1, ncol=3,labels = c("A", "B", "C"), label_size = 18)
    title <- ggdraw() + draw_label(paste0(fig_title), size = 10)
    plot_grid(title, pic, legend, ncol=1, rel_heights=c(0.3, 1.2, 0.6),scale = c(1, 1,1)) # rel_heights values control title margins

}

data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))
pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_spe(x[1]) )
pic <- cowplot::plot_grid(pic_all[[1]], pic_all[[2]], pic_all[[3]],nrow = 3, ncol=1)

ggsave(pic, filename = file.path(figure_dir, "dm6.performance_by_depth.pdf"), width = 7.5, height = 8, bg = 'transparent')





performance_plot_spe <- function(te_class){

    genome <- "dm6"
    data_path <- file.path(result_root, "simulation", genome)


    # perfomance <- read.table(paste0(data_path,"/intersection/pp_intersection/data/performance.txt"))
    perfomance <- read.table(paste0(data_path,"/intersection/pp_intersection/data/performance.",te_class,".txt"))
    # print(ncol(perfomance))
    
    if(ncol(perfomance) != 6){
        return(1)
    }
    colnames(perfomance) <- c("data_type", "depth", "caller", "Sensitivity", "Precision", "F1-score")
    perfomance <- perfomance[ perfomance$caller!="TLDR_all" & perfomance$caller!="LOCATE_all", ]

    caller_dm6_dis <- c("TEMP2", "MELT","TLDR_pass","TELR", "TrEMOLO", "GraffiTE","LOCATE_pass")

    perfomance$caller <- factor(perfomance$caller, levels=caller_dm6_dis)
    perfomance <- perfomance[ ! ( perfomance$caller=="TELR" & perfomance$Precision=='0' ), ]
    perfomance <- perfomance[ ! ( perfomance$caller=="MELT" & perfomance$Precision=='0' ), ]
    perfomance <- perfomance[ ! ( perfomance$caller=="TrEMOLO" & perfomance$Precision=='0' ), ]

    perfomance <- perfomance[ perfomance$caller!="TLDR_all" & perfomance$caller!="xTea_all" & perfomance$caller!="LOCATE_all", ]




    data_type <- "simulation_germ_clr"

    pre_data <- perfomance[perfomance$data_type==data_type,]       

    p1 <- performance_plot(pre_data, "Sensitivity","None")
    p2 <- performance_plot(pre_data, "Precision","None")
    p3 <- performance_plot(pre_data, "F1-score","None")
    p_legend <- performance_plot(pre_data, "Sensitivity","top")
    legend <- get_plot_component(p_legend, 'guide-box-top', return_all = TRUE)

    pic <- cowplot::plot_grid(p1, p2, p3,nrow = 1, ncol=3)
    title <- ggdraw() + draw_label(paste0(te_class), fontface='bold')
    plot_grid(title, pic, legend, ncol=1, rel_heights=c(0.15, 1, 0.3),scale = c(1, 1,1)) # rel_heights values control title margins

}

# data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))
data_type_matrix <- data.frame(TE_CLASS=c('LINE', 'DNA', 'LTR', "RC")) #

pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_spe(x[1]) )
pic <- cowplot::plot_grid(pic_all[[1]], pic_all[[2]],  pic_all[[3]], pic_all[[4]], nrow = 4, ncol=1, labels = c("A", "B", "C", "D")) #

ggsave(pic, filename = file.path(figure_dir, "dm6.performance_by_TE_class.pdf"), width = 6, height = 8)


### distance divergency
genome <- "dm6"
data_path <- file.path(result_root, "simulation", genome)

perfomance_dis <- read.table(paste0(data_path, "/intersection/all.merge.dis.txt"), header = T)
# perfomance_dis <- perfomance_dis[perfomance_dis$frequency>0.2,]
head(perfomance_dis)
# bp_dis
perfomance.bp_dis <- melt(perfomance_dis, id.vars=c("insID", "sample", "depth"), measure.vars =c("LOCATE_all.bp_dis", "LOCATE_pass.bp_dis", "MELT.bp_dis","TEMP2.bp_dis", "TLDR_all.bp_dis","TLDR_pass.bp_dis" , "TELR.bp_dis","TrEMOLO.bp_dis","GraffiTE.bp_dis"), variable.name = "caller", value.name = "bp_dis")
perfomance.bp_dis$caller <- sub("\\..*", "", perfomance.bp_dis$caller)
perfomance.bp_dis <- perfomance.bp_dis[perfomance.bp_dis$bp_dis != "-", ]
perfomance.bp_dis$bp_dis <- as.numeric(perfomance.bp_dis$bp_dis)
# te_dis
perfomance.te_dis <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.te_dis", "LOCATE_pass.te_dis", "MELT.te_dis","TEMP2.te_dis", "TLDR_all.te_dis","TLDR_pass.te_dis" , "TELR.te_dis","TrEMOLO.te_dis","GraffiTE.te_dis" ), variable.name = "caller", value.name = "te_dis")
perfomance.te_dis$caller <- sub("\\..*", "", perfomance.te_dis$caller)
perfomance.te_dis <- perfomance.te_dis[perfomance.te_dis$te_dis != "-", ]
perfomance.te_dis$te_dis <- as.numeric(perfomance.te_dis$te_dis)

# te_identity
perfomance.te_identity <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.te_identity", "LOCATE_pass.te_identity", "MELT.te_identity","TEMP2.te_identity", "TLDR_all.te_identity","TLDR_pass.te_identity" , "TELR.te_identity","TrEMOLO.te_identity","GraffiTE.te_identity" ), variable.name = "caller", value.name = "te_identity")
perfomance.te_identity$caller <- sub("\\..*", "", perfomance.te_identity$caller)
perfomance.te_identity <- perfomance.te_identity[perfomance.te_identity$te_identity != "-", ]
perfomance.te_identity$te_identity <- as.numeric(perfomance.te_identity$te_identity)


# div
perfomance.te_len <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.te_len", "LOCATE_pass.te_len", "MELT.te_len","TEMP2.te_len", "TLDR_all.te_len","TLDR_pass.te_len" , "TELR.te_len","TrEMOLO.te_len","GraffiTE.te_len" ), variable.name = "caller", value.name = "te_len")
perfomance.te_len$caller <- sub("\\..*", "", perfomance.te_len$caller)
perfomance.div <- melt(perfomance_dis, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_all.div", "LOCATE_pass.div", "MELT.div","TEMP2.div", "TLDR_all.div","TLDR_pass.div" , "TELR.div","TrEMOLO.div","GraffiTE.div" ), variable.name = "caller", value.name = "div")
perfomance.div$caller <- sub("\\..*", "", perfomance.div$caller)
perfomance.len_div <- merge(perfomance.te_len, perfomance.div, by = c("insID","sample", "depth", "caller"))
perfomance.len_div <- perfomance.len_div[perfomance.len_div$div != "-", ]
perfomance.len_div$div <- as.numeric(perfomance.len_div$div)
perfomance.len_div$te_len <- as.numeric(perfomance.len_div$te_len)

# LOCATE_TLDR
perfomance_dis.LOCATE_TLDR <- perfomance_dis[perfomance_dis$LOCATE_pass.div!="-" & perfomance_dis$TLDR_pass.div!="-", ]
perfomance.LOCATE_TLDR.te_len <- melt(perfomance_dis.LOCATE_TLDR, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_pass.te_len","TLDR_pass.te_len"  ), variable.name = "caller", value.name = "te_len")
perfomance.LOCATE_TLDR.te_len$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR.te_len$caller)
perfomance.LOCATE_TLDR.div <- melt(perfomance_dis.LOCATE_TLDR, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_pass.div","TLDR_pass.div" ), variable.name = "caller", value.name = "div")
perfomance.LOCATE_TLDR.div$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR.div$caller)
perfomance.LOCATE_TLDR.len_div <- merge(perfomance.LOCATE_TLDR.te_len, perfomance.LOCATE_TLDR.div, by = c("insID","sample", "depth", "caller"))
perfomance.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[perfomance.LOCATE_TLDR.len_div$div != "-", ]
perfomance.LOCATE_TLDR.len_div$div <- as.numeric(perfomance.LOCATE_TLDR.len_div$div)
perfomance.LOCATE_TLDR.len_div$te_len <- as.numeric(perfomance.LOCATE_TLDR.len_div$te_len)

# LOCATE_TLDR
perfomance_dis.LOCATE_TLDR_TELR <- perfomance_dis[perfomance_dis$LOCATE_pass.div!="-" & perfomance_dis$TLDR_pass.div!="-" & perfomance_dis$TELR.div!="-", , ]
perfomance.LOCATE_TLDR_TELR.te_len <- melt(perfomance_dis.LOCATE_TLDR_TELR, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_pass.te_len","TLDR_pass.te_len" ,"TELR.te_len"  ), variable.name = "caller", value.name = "te_len")
perfomance.LOCATE_TLDR_TELR.te_len$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_TELR.te_len$caller)
perfomance.LOCATE_TLDR_TELR.div <- melt(perfomance_dis.LOCATE_TLDR_TELR, id.vars=c("insID","sample", "depth"), measure.vars =c("LOCATE_pass.div","TLDR_pass.div","TELR.div" ), variable.name = "caller", value.name = "div")
perfomance.LOCATE_TLDR_TELR.div$caller <- sub("\\..*", "", perfomance.LOCATE_TLDR_TELR.div$caller)
perfomance.LOCATE_TLDR_TELR.len_div <- merge(perfomance.LOCATE_TLDR_TELR.te_len, perfomance.LOCATE_TLDR_TELR.div, by = c("insID","sample", "depth", "caller"))
perfomance.LOCATE_TLDR_TELR.len_div <- perfomance.LOCATE_TLDR_TELR.len_div[perfomance.LOCATE_TLDR_TELR.len_div$div != "-", ]
perfomance.LOCATE_TLDR_TELR.len_div$div <- as.numeric(perfomance.LOCATE_TLDR_TELR.len_div$div)
perfomance.LOCATE_TLDR_TELR.len_div$te_len <- as.numeric(perfomance.LOCATE_TLDR_TELR.len_div$te_len)


perfomance.bp_dis <- perfomance.bp_dis[! perfomance.bp_dis$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.te_dis <- perfomance.te_dis[! perfomance.te_dis$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.te_identity <- perfomance.te_identity[! perfomance.te_identity$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.len_div <- perfomance.len_div[! perfomance.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[! perfomance.LOCATE_TLDR.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]
perfomance.LOCATE_TLDR_TELR.len_div <- perfomance.LOCATE_TLDR_TELR.len_div[! perfomance.LOCATE_TLDR_TELR.len_div$caller %in% c("TLDR_all",  "xTea_all", "LOCATE_all"),]


# levels
caller_dm6_dis <- c("TEMP2","MELT", "TLDR_pass","TELR","TrEMOLO","GraffiTE","LOCATE_pass")
perfomance.bp_dis$caller <- factor(perfomance.bp_dis$caller, levels=caller_dm6_dis)
perfomance.te_dis$caller <- factor(perfomance.te_dis$caller, levels=caller_dm6_dis)
perfomance.te_identity$caller <- factor(perfomance.te_identity$caller, levels=caller_dm6_dis)
perfomance.len_div$caller <- factor(perfomance.len_div$caller, levels=caller_dm6_dis)
perfomance.LOCATE_TLDR.len_div$caller <- factor(perfomance.LOCATE_TLDR.len_div$caller, levels=caller_dm6_dis)
perfomance.LOCATE_TLDR_TELR.len_div$caller <- factor(perfomance.LOCATE_TLDR_TELR.len_div$caller, levels=caller_dm6_dis)



# plot divergency

performance_plot_dis_spe <- function(data_type){
    pre_data.bp_dis <- perfomance.bp_dis[perfomance.bp_dis$sample==data_type,]
    pre_data.bp_dis <- pre_data.bp_dis %>% group_by(caller, depth) %>% dplyr::summarise(avg_bp_dis = mean(bp_dis), .groups = "keep" )

    pre_data.te_identity <- perfomance.te_identity[perfomance.te_identity$sample==data_type,]
    pre_data.te_identity <- pre_data.te_identity %>% group_by(caller, depth) %>% dplyr::summarise(avg_te_identity = mean(te_identity), .groups = "keep")

    pre_data.len_div <- perfomance.len_div[perfomance.len_div$sample==data_type,]

    pre_data.len_div <- pre_data.len_div %>% group_by(caller, depth) %>% dplyr::summarise(avg_div = sum(te_len*div)/sum(te_len), .groups = "keep")
    

    pre_data.len_div$avg_div <- log(pre_data.len_div$avg_div, 10)

    pre_data.LOCATE_TLDR.len_div <- perfomance.LOCATE_TLDR.len_div[perfomance.LOCATE_TLDR.len_div$sample==data_type,]
    pre_data.LOCATE_TLDR.len_div <- pre_data.LOCATE_TLDR.len_div %>% group_by(caller, depth) %>% dplyr::summarise(avg_div = sum(te_len*div)/sum(te_len), .groups = "keep")
    pre_data.LOCATE_TLDR.len_div.1 <-  dcast( pre_data.LOCATE_TLDR.len_div, depth~caller, value.var = "avg_div" )

    pre_data.LOCATE_TLDR.len_div$avg_div <- log(pre_data.LOCATE_TLDR.len_div$avg_div, 10)

    
    # write.table(pre_data.LOCATE_TLDR_xTea.len_div, paste0(data_path, "/figure/",genome,".",data_type, ".div.txt"),row.names = FALSE, sep ="\t",col.names =FALSE, quote =FALSE)

    # p3 <- performance_plot_dis(pre_data.len_div, "avg_div","None",data_type,"All insertion" , "FALSE")
    # p4 <- performance_plot_dis(pre_data.LOCATE_TLDR.len_div, "avg_div","None",data_type)
    p5 <- performance_plot_dis(pre_data.LOCATE_TLDR.len_div, "avg_div","None",data_type,"All insertion" ,"FALSE")

    p_legend <- performance_plot_dis(pre_data.LOCATE_TLDR.len_div, "avg_div","top",data_type)
    legend <- get_plot_component(p_legend, 'guide-box-top', return_all = TRUE)
    

    # pic <- cowplot::plot_grid(p3, p5, nrow = 1, ncol=2, rel_widths = c(1.15,1))
    # pic <- cowplot::plot_grid(p5, nrow = 1, ncol=1)

    if (data_type == "simulation_germ_ccs"){
        fig_title <- "HiFi"
    }
    if (data_type == "simulation_germ_clr"){
        fig_title <- "CLR"
    }
    if (data_type == "simulation_germ_ont"){
        fig_title <- "ONT"
    }

    title <- ggdraw() + draw_label(paste0(fig_title), size=10)
    p.tmp <- plot_grid(title, p5, ncol=1, rel_heights=c(0.3, 1.2)) # rel_heights values control title margins
    return(list(p.tmp, legend, pre_data.LOCATE_TLDR.len_div.1 ))
}

data_type_matrix <- data.frame(data_type=c('simulation_germ_ccs', 'simulation_germ_clr', 'simulation_germ_ont'))
pic_all <- apply(data_type_matrix, 1, function(x) performance_plot_dis_spe(x[1]) )
pic <- cowplot::plot_grid(pic_all[[1]][[1]], pic_all[[2]][[1]],pic_all[[3]][[1]],nrow = 1, ncol=3)
# pic.1 <- cowplot::plot_grid(pic, pic_all[[1]][[2]],nrow = 2, ncol=1, rel_heights=c(1, 0.4),scale = c(1,1))

ggsave(pic, filename = file.path(figure_dir, "dm6.sequence_divergence_by_depth.pdf"), width = 6, height = 2.4)

len_div.hifi <-  pic_all[[1]][[3]]
len_div.hifi$data_type <- "simulation_germ_ccs"
len_div.clr <-  pic_all[[2]][[3]]
len_div.clr$data_type <- "simulation_germ_clr"
len_div.ont <-  pic_all[[3]][[3]]
len_div.ont$data_type <- "simulation_germ_ont"

len_div <- rbind(len_div.hifi, len_div.clr, len_div.ont)

write.xlsx(len_div, file.path(table_dir, "dm6.sequence_divergence_summary.xlsx"), sheetName="dm6.dis", colNames = TRUE, append = TRUE)

            


