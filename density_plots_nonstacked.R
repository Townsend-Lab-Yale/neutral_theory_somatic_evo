
# function to create different density plots ---- 
density_plotter_function <- function(main_df,
                                     this_tumor,
                                     main_size=10,
                                     sub_1 = "PIK3CA E545K",
                                     sub_1_lab_y = 0.1,
                                     sub_1_lab_x = 1,
                                     sub_1_arrow_x = 1,
                                     sub_1_arrow_y = 0.01,
                                     sub_1_lab_lab = sub_1,
                                     sub_1_arrow_alpha=1,
                                     sub_2,
                                     sub_2_lab_y = 0.1,
                                     sub_2_lab_x = 1,
                                     sub_2_arrow_x = 1,
                                     sub_2_arrow_y = 0.01,
                                     sub_2_lab_lab = sub_2,
                                     sub_2_arrow_alpha = 1,
                                     point_shape = "|",
                                     point_size = 10,
                                     x_axis_breaks = c(1,1e4,1e5),
                                     x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5)),
                                     y_axis_breaks = c(0,.1,.2),
                                     y_axis_labels = c(0,.1,.2),
                                     y_axis_limits = c(0,.25),
                                     figure_title=this_tumor,
                                     extender_percent=0.01,
                                     minimum_maximum=1.1e6,
                                     margins=c(2,2,2,2)
){
  # recurrent mutations in this tumor
  data_for_this_tumor_recur <- subset(main_df, tumor_type==this_tumor & freq>1)
  
  # data for substitution 1 
  sub_1_data <- subset(data_for_this_tumor_recur, name==sub_1)
  message(paste("Selection intensity of substitution 1:",round(sub_1_data$gamma_epistasis,2)))
  sub_1_df <- data.frame(sub_1_lab_x,
                         sub_1_lab_y,
                         sub_1_arrow_x=sub_1_arrow_x,
                         sub_1_arrow_y,
                         sub_1_lab_lab)
  
  # data for substitution 2 
  sub_2_data <- subset(data_for_this_tumor_recur, name==sub_2)
  message(paste("Selection intensity of substitution 2:",round(sub_2_data$gamma_epistasis,2)))
  sub_2_df <- data.frame(sub_2_lab_x,
                         sub_2_lab_y,
                         sub_2_arrow_x=sub_2_arrow_x,
                         sub_2_arrow_y,
                         sub_2_lab_lab)
  
  
  # create the plot 
  this_density_plot <- ggplot(data = subset(main_df, tumor_type==this_tumor), 
                              aes(x=gamma_epistasis)) + 
    geom_density(aes(fill=synonymous),alpha=.5) + 
    
    # first substitution 
    geom_point(data = sub_1_data, 
               aes(x=gamma_epistasis,
                   y=c(0)),
               shape=point_shape,size=point_size) + 
    geom_text(data=sub_1_df,
              aes(x=sub_1_lab_x,
                  y=sub_1_lab_y),
              label=sub_1_lab_lab,
              size=main_size*(5/14),
              vjust=0) +
    geom_segment(data=sub_1_df,
                 aes(x=sub_1_lab_x,
                     y=sub_1_lab_y,
                     xend=sub_1_arrow_x,
                     yend=sub_1_arrow_y),
                 arrow = arrow(length = unit(0.05, "npc")),
                 alpha=sub_1_arrow_alpha) + 
    
    
    # second substitution 
    geom_point(data = sub_2_data, 
               aes(x=gamma_epistasis,
                   y=c(0)),
               shape=point_shape,size=point_size) + 
    geom_text(data=sub_2_df,
              aes(x=sub_2_lab_x,
                  y=sub_2_lab_y),
              label=sub_2_lab_lab,
              size=main_size*(5/14),
              vjust=0) +
    geom_segment(data=sub_2_df,
                 aes(x=sub_2_lab_x,
                     y=sub_2_lab_y,
                     xend=sub_2_arrow_x,
                     yend=sub_2_arrow_y),
                 arrow = arrow(length = unit(0.05, "npc")),
                 alpha=sub_2_arrow_alpha) + 
    # title in the plot
    
    # annotate("text",label=this_tumor,x=(max(data_for_this_tumor_recur$gamma_epistasis)+(max(data_for_this_tumor_recur$gamma_epistasis)*extender_percent))/2,y=)
    
    # scaling, etc. 
    scale_x_continuous(breaks=x_axis_breaks,
                       labels=x_axis_labels,
                       trans="mysqrt",
                       expand = c(0,0)) + 
    coord_cartesian(xlim = c(1,ifelse(minimum_maximum>max(data_for_this_tumor_recur$gamma_epistasis),minimum_maximum+(minimum_maximum*extender_percent),
                                      (max(data_for_this_tumor_recur$gamma_epistasis)+(max(data_for_this_tumor_recur$gamma_epistasis)*extender_percent))))) + 
    guides(fill=FALSE)+
    scale_y_continuous(expand=c(0,0), 
                       limits=y_axis_limits,
                       breaks=y_axis_breaks,
                       labels = y_axis_labels) + 
    theme(  axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_text(size=main_size),
            text = element_text(size = main_size),
            plot.title = element_text(size=main_size,face = "plain"),
            plot.margin=unit(c(0,0,0,0), "mm")) #+
  # theme_classic() + 
  # labs(title=figure_title) 
  
  
  return(this_density_plot)
}



# test_plot <- density_plotter_function(main_df = combined_all_data.noNA,this_tumor = "LUAD",sub_1 = "KRAS G12C")
# 
# 
# test_plot <- density_plotter_function(main_df = combined_all_data.noNA,this_tumor = "HNSC_HPVneg",sub_1 = "PIK3CA E545K",figure_title = "HPV negative HNSCC",sub_1_lab_x = 5e4,sub_2 = "TP53 R283P",sub_2_lab_x = 1e5,sub_2_lab_y = .2,main_size = 10,sub_2_arrow_y = .04,sub_1_arrow_y = 0.04,extender_percent = 0.02,sub_2_arrow_x = 358905.26-40000,sub_1_arrow_x = 13434.61-1000)
# 
# 
# library(cowplot)
# 
# test_plot_combine <- plot_grid(test_plot,test_plot,test_plot,test_plot,test_plot,test_plot,ncol=3,align = 'hv')
# 
# save_plot(filename = "figures/test_combined.png",base_width = 6.5,plot = test_plot_combine,dpi=600)






# building plots 

# common parameters 

common.size <-7

# BLCA ----
BLCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "BLCA",
                                         main_size = common.size,
                                         sub_1 = "FBXW7 R425G",
                                         sub_2 = "HRAS G13R",
                                         sub_1_lab_y = .12,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x = 372257.26,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 7e4,
                                         sub_2_arrow_x = 96564.92,
                                         sub_2_arrow_y = .03,
                                         extender_percent = 0.06,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

BLCA_density <- ggdraw(BLCA_density) + draw_label("BLCA",x = .5,y=.85,size = common.size)

save_plot(filename = "figures/BLCA_density.png",plot = BLCA_density,base_height = 5/6,base_width = 6.5/4)

# BRCA ----
BRCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "BRCA",
                                         main_size = common.size,
                                         sub_1 = "PIK3CA H1047R",
                                         sub_2 = "TP53 Y220S",
                                         sub_1_lab_y = .05,
                                         sub_1_lab_x = 6e5,
                                         sub_1_arrow_x =226324.32,
                                         sub_1_arrow_y = 0.03,
                                         sub_2_lab_y = .12,
                                         sub_2_lab_x = 7e4,
                                         sub_2_arrow_x = 81082.88,
                                         sub_2_arrow_y = .03,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

BRCA_density <- ggdraw(BRCA_density) + draw_label("BRCA",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/BRCA_density.png",plot = BRCA_density,base_height = 5/6,base_width = 6.5/4)


# CESC ---- 

CESC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "CESC",
                                         main_size = common.size,
                                         sub_1 = "FBXW7 R505G",
                                         sub_2 = "MED12 D23Y",
                                         sub_1_lab_y = .1,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x =531100.12,
                                         sub_1_arrow_y = 0.03,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x = 68681.85,
                                         sub_2_arrow_y = .03,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

CESC_density <- ggdraw(CESC_density) + draw_label("CESC",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/CESC_density.png",plot = CESC_density,base_height = 5/6,base_width = 6.5/4)


# COAD ----

COAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "COAD",
                                         main_size = common.size,
                                         sub_1 = "BRAF V600E",
                                         sub_2 = "TP53 A159P",
                                         sub_1_lab_y = .2,
                                         sub_1_lab_x = 4e5,
                                         sub_1_arrow_x =1200000.5,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .1,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x = 369823.27,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         extender_percent = 0.02,
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4)

COAD_density <- ggdraw(COAD_density) + draw_label("COAD",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/COAD_density.png",plot = COAD_density,base_height = 5/6,base_width = 6.5/4)

# ESCA ---- 


ESCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "ESCA",
                                         main_size = common.size,
                                         sub_1 = "TP53 R175G",
                                         sub_2 = "PIK3CA H1047L",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 7e5,
                                         sub_1_arrow_x = 1707238.91,
                                         sub_1_arrow_y = 0.02,
                                         sub_2_lab_y = .045,
                                         sub_2_lab_x = 240937.09,
                                         sub_2_arrow_x = 240937.09,
                                         sub_2_arrow_y = .02,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                         y_axis_limits = c(0,.14),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

ESCA_density <- ggdraw(ESCA_density) + draw_label("ESCA",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/ESCA_density.png",plot = ESCA_density,base_height = 5/6,base_width = 6.5/4)


# GBM ----

GBM_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                        this_tumor = "GBM",
                                        main_size = common.size,
                                        sub_1 = "BRAF V600E",
                                        sub_2 = "TP53 Y220C",
                                        sub_1_lab_y = .04,
                                        sub_1_lab_x = 5e5,
                                        sub_1_arrow_x = 60377.811,
                                        sub_1_arrow_y = 0.015,
                                        sub_2_lab_y = .08,
                                        sub_2_lab_x = 5e4,
                                        sub_2_arrow_x =  32406.12,
                                        sub_2_arrow_y = .03,
                                        x_axis_breaks = c(1,1e4,1e5,1e6),x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                        y_axis_limits = c(0,.15),
                                        y_axis_breaks = c(0,0.1),
                                        y_axis_labels = c(0,0.1),
                                        point_size = 4)

GBM_density <- ggdraw(GBM_density) + draw_label("GBM",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/GBM_density.png",plot = GBM_density,base_height = 5/6,base_width = 6.5/4)



# HNSC_HPVpos ----

HNSC_HPVpos_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                                this_tumor = "HNSC_HPVpos",
                                                main_size = common.size,
                                                sub_1 = "FBXW7 R505G",
                                                sub_2 = "FGFR3 S249C",
                                                sub_1_lab_y = .06,
                                                sub_1_lab_x = 4e5,
                                                sub_1_arrow_x = 1140442.02,
                                                sub_1_arrow_y = 0.02,
                                                sub_2_lab_y = .03,
                                                sub_2_lab_x = 2e5,
                                                sub_2_arrow_x =  145184.52,
                                                sub_2_arrow_y = .02,
                                                x_axis_breaks = c(1,1e4,1e5,1e6),
                                                x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                                sub_1_arrow_alpha = 1,
                                                sub_2_arrow_alpha = 1,
                                                y_axis_limits = c(0,.11),
                                                y_axis_breaks = c(0,0.1),
                                                y_axis_labels = c(0,0.1),
                                                point_size = 4)

HNSC_HPVpos_density <- ggdraw(HNSC_HPVpos_density) + draw_label(label = expression(HPV^{"+"}~HNSC),x = .55,y=.9,size = common.size)

save_plot(filename = "figures/HNSC_HPVpos_density.png",plot = HNSC_HPVpos_density,base_height = 5/6,base_width = 6.5/4)


# HNSC_HPVneg ----

HNSC_HPVneg_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                                this_tumor = "HNSC_HPVneg",
                                                main_size = common.size,
                                                sub_1 = "TP53 R283P",
                                                sub_2 = "TP53 E298*",
                                                sub_1_lab_y = .12,
                                                sub_1_lab_x = 5e5,
                                                sub_1_arrow_x = 358905.26,
                                                sub_1_arrow_y = 0.05,
                                                sub_2_lab_y = .08,
                                                sub_2_lab_x = 6e4,
                                                sub_2_arrow_x =  200186.43,
                                                sub_2_arrow_y = .04,
                                                x_axis_breaks = c(1,1e4,1e5,1e6),
                                                x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                                extender_percent = 0.06,
                                                sub_1_arrow_alpha = 1,
                                                sub_2_arrow_alpha = 1,
                                                y_axis_limits = c(0,.22),
                                                y_axis_breaks = c(0,0.1,0.2),
                                                y_axis_labels = c(0,0.1,0.2),
                                                point_size = 4)

HNSC_HPVneg_density <- ggdraw(HNSC_HPVneg_density) + draw_label(label = expression(HPV^{"−"}~HNSC),x = .55,y=.9,size = common.size)

save_plot(filename = "figures/HNSC_HPVneg_density.png",plot = HNSC_HPVneg_density,base_height = 5/6,base_width = 6.5/4)


# KIRC ---- 


KIRC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "KIRC",
                                         main_size = common.size,
                                         sub_1 = "VHL S68*",
                                         sub_2 = "VHL S65*",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  185360.85,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .05,
                                         sub_2_lab_x = 8.5e4,
                                         sub_2_arrow_x =  118487.13,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.14),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

KIRC_density <- ggdraw(KIRC_density) + draw_label(label = "KIRC",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/KIRC_density.png",plot = KIRC_density,base_height = 5/6,base_width = 6.5/4)


# LAML ----

LAML_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LAML",
                                         main_size = common.size,
                                         sub_1 = "KIT D816V",
                                         sub_2 = "KRAS G12V",
                                         sub_1_lab_y = .065,
                                         sub_1_lab_x = 8e5,
                                         sub_1_arrow_x =  1969488.61,
                                         sub_1_arrow_y = 0.02,
                                         sub_2_lab_y = .04,
                                         sub_2_lab_x = 3e5,
                                         sub_2_arrow_x =   474810.3,
                                         sub_2_arrow_y = .02,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.11),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

LAML_density <- ggdraw(LAML_density) + draw_label(label = "LAML",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LAML_density.png",plot = LAML_density,base_height = 5/6,base_width = 6.5/4)


# LGG ---- 
LGG_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                        this_tumor = "LGG",
                                        main_size = common.size,
                                        sub_1 = "IDH1 R132G",
                                        sub_2 = "IDH1 R132S",
                                        sub_1_lab_y = .08,
                                        sub_1_lab_x = 6e6,
                                        sub_1_arrow_x =  13233018.63,
                                        sub_1_arrow_y = 0.025,
                                        sub_2_lab_y = .04,
                                        sub_2_lab_x = 1e6,
                                        sub_2_arrow_x =   2848965.8,
                                        sub_2_arrow_y = .025,
                                        x_axis_breaks = c(1,1e5,1e6,1e7),
                                        x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.02,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                        y_axis_limits = c(0,.15),
                                        y_axis_breaks = c(0,0.1),
                                        y_axis_labels = c(0,0.1),
                                        point_size = 4)

LGG_density <- ggdraw(LGG_density) + draw_label(label = "LGG",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LGG_density.png",plot = LGG_density,base_height = 5/6,base_width = 6.5/4)



# LIHC ---- 
LIHC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LIHC",
                                         main_size = common.size,
                                         sub_1 = "BAZ2A R1642L",
                                         sub_2 = "CTNNB1 I35S",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 3e5,
                                         sub_1_arrow_x =  949298.36,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .04,
                                         sub_2_lab_x = 1.5e5,
                                         sub_2_arrow_x =   228601.05,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1
                                         ,sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.15),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

LIHC_density <- ggdraw(LIHC_density) + draw_label(label = "LIHC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LIHC_density.png",plot = LIHC_density,base_height = 5/6,base_width = 6.5/4)

# LUAD ---- 

LUAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LUAD",
                                         main_size = common.size,
                                         sub_1 = "EGFR L858R",
                                         sub_2 = "CTNNB1 S37F",
                                         sub_1_lab_y = .15,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  280236.22,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .08,
                                         sub_2_lab_x = 6e4,
                                         sub_2_arrow_x =    120122.92,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LUAD_density <- ggdraw(LUAD_density) + draw_label(label = "LUAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LUAD_density.png",plot = LUAD_density,base_height = 5/6,base_width = 6.5/4)




# LUSC ---- 

LUSC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "LUSC",
                                         main_size = common.size,
                                         sub_1 = "TP53 Y234S",
                                         sub_2 = "NFE2L2 R34G",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =  90448.39,
                                         sub_1_arrow_y = 0.03,
                                         sub_2_lab_y = .13,
                                         sub_2_lab_x =5e4,
                                         sub_2_arrow_x =    69242.49,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.22),
                                         y_axis_breaks = c(0,0.1,0.2),
                                         y_axis_labels = c(0,0.1,0.2),
                                         point_size = 4)

LUSC_density <- ggdraw(LUSC_density) + draw_label(label = "LUSC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/LUSC_density.png",plot = LUSC_density,base_height = 5/6,base_width = 6.5/4)




# OV ---- 

OV_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                       this_tumor = "OV",
                                       main_size = common.size,
                                       sub_1 = "MGRN1 G64A",
                                       sub_2 = "TP53 V157F",
                                       sub_1_lab_y = .12,
                                       sub_1_lab_x = 5e5,
                                       sub_1_arrow_x =   409080.25,
                                       sub_1_arrow_y = 0.025,
                                       sub_2_lab_y = .08,
                                       sub_2_lab_x = 8e4,
                                       sub_2_arrow_x =    184993.66,
                                       sub_2_arrow_y = .025,
                                       x_axis_breaks = c(1,1e4,1e5,1e6),
                                       x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                       sub_1_arrow_alpha = 1,
                                       sub_2_arrow_alpha = 1,
                                       y_axis_limits = c(0,.22),
                                       y_axis_breaks = c(0,0.1,0.2),
                                       y_axis_labels = c(0,0.1,0.2),
                                       point_size = 4)

OV_density <- ggdraw(OV_density) + draw_label(label = "OV",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/OV_density.png",plot = OV_density,base_height = 5/6,base_width = 6.5/4)



# PAAD ---- 

PAAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "PAAD",
                                         main_size = common.size,
                                         sub_1 = "KRAS G12R",
                                         sub_2 = "KRAS G12V",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 2e6,
                                         sub_1_arrow_x =   9513627.87,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .04,
                                         sub_2_lab_x = 8e5,
                                         sub_2_arrow_x =    2801115.02,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6)),extender_percent = 0.02,sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.14),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

PAAD_density <- ggdraw(PAAD_density) + draw_label(label = "PAAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/PAAD_density.png",plot = PAAD_density,base_height = 5/6,base_width = 6.5/4)


# PRAD ---- 

PRAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "PRAD",
                                         main_size = common.size,
                                         sub_1 = "SPOP Y87S",
                                         sub_2 = "SPOP W131G",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =   683195.22,
                                         sub_1_arrow_y = 0.025,
                                         sub_2_lab_y = .04,
                                         sub_2_lab_x = 1e5,
                                         sub_2_arrow_x =    552894.51,
                                         sub_2_arrow_y = .025,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.06,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.14),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

PRAD_density <- ggdraw(PRAD_density) + draw_label(label = "PRAD",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/PRAD_density.png",plot = PRAD_density,base_height = 5/6,base_width = 6.5/4)


# READ ---- 

READ_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "READ",
                                         main_size = common.size,
                                         sub_1 = "APC E1322*",
                                         sub_2 = "APC E1209*",
                                         sub_1_lab_y = .08,
                                         sub_1_lab_x = 2e5,
                                         sub_1_arrow_x =   2613352.63,
                                         sub_1_arrow_y = 0.02,
                                         sub_2_lab_y = .03,
                                         sub_2_lab_x = 2e5,
                                         sub_2_arrow_x =    2513352.63,
                                         sub_2_arrow_y = .01,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),extender_percent = 0.05,sub_1_arrow_alpha = 1,sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.15),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

READ_density <- ggdraw(READ_density) + draw_label(label = "READ",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/READ_density.png",plot = READ_density,base_height = 5/6,base_width = 6.5/4)



# SKCMM ---- 

SKCMM_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                          this_tumor = "SKCMM",
                                          main_size = common.size,
                                          sub_1 = "NRAS Q61R",
                                          sub_2 = "NRAS Q61K",
                                          sub_1_lab_y = .2,
                                          sub_1_lab_x = 2e6,
                                          sub_1_arrow_x =   4707860.13,
                                          sub_1_arrow_y = 0.05,
                                          sub_2_lab_y = .1,
                                          sub_2_lab_x = 1e6,
                                          sub_2_arrow_x =    4438279.31,
                                          sub_2_arrow_y = .02,
                                          x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                          x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                          extender_percent = 0.02,
                                          sub_1_arrow_alpha = 1,
                                          sub_2_arrow_alpha = 1,
                                          y_axis_limits = c(0,.4),
                                          y_axis_breaks = c(0,0.1,0.2,0.3),
                                          y_axis_labels = c(0,0.1,0.2,0.3),
                                          point_size = 4)

SKCMM_density <- ggdraw(SKCMM_density) + draw_label(label = "SKCMM",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/SKCMM_density.png",plot = SKCMM_density,base_height = 5/6,base_width = 6.5/4)



# SKCMP ---- 

SKCMP_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                          this_tumor = "SKCMP",
                                          main_size = common.size,
                                          sub_1 = "BRAF V600E",
                                          sub_2 = "KIT K642E",
                                          sub_1_lab_y = .1,
                                          sub_1_lab_x = 1e6,
                                          sub_1_arrow_x =   3092084.99,
                                          sub_1_arrow_y = 0.03,
                                          sub_2_lab_y = .05,
                                          sub_2_lab_x = 1e5,
                                          sub_2_arrow_x =    127735.46,
                                          sub_2_arrow_y = .03,
                                          x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                          x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                          extender_percent = 0.02,
                                          sub_1_arrow_alpha = 1,
                                          sub_2_arrow_alpha = 1,
                                          y_axis_limits = c(0,.22),
                                          y_axis_breaks = c(0,0.1,0.2),
                                          y_axis_labels = c(0,0.1,0.2),
                                          point_size = 4)

SKCMP_density <- ggdraw(SKCMP_density) + draw_label(label = "SKCMP",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/SKCMP_density.png",plot = SKCMP_density,base_height = 5/6,base_width = 6.5/4)



# STAD ----

STAD_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "STAD",
                                         main_size = common.size,
                                         sub_1 = "RHOA Y42S",
                                         sub_2 = "PIK3CA H1047R",
                                         sub_1_lab_y = .13,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =65728.42,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .25,
                                         sub_2_lab_x = 5e4,
                                         sub_2_arrow_x = 27539.21,
                                         sub_2_arrow_y = .06,
                                         x_axis_breaks = c(1,1e4,1e5,1e6),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6)),
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4,extender_percent = 0.06)

STAD_density <- ggdraw(STAD_density) + draw_label("STAD",x = .5,y=.9,size = common.size)

save_plot(filename = "figures/STAD_density.png",plot = STAD_density,base_height = 5/6,base_width = 6.5/4)


# THCA ---- 

THCA_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "THCA",
                                         main_size = common.size,
                                         sub_1 = "BRAF V600E",
                                         sub_2 = "NRAS Q61R",
                                         sub_1_lab_y = .06,
                                         sub_1_lab_x = 1e7,
                                         sub_1_arrow_x =   33029682.79,
                                         sub_1_arrow_y = 0.02,
                                         sub_2_lab_y = .04,
                                         sub_2_lab_x = 2e6,
                                         sub_2_arrow_x =    4601059.2,
                                         sub_2_arrow_y = .02,
                                         x_axis_breaks = c(1,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^5),expression(10^6),expression(10^7)),extender_percent = 0.02,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.11),
                                         y_axis_breaks = c(0,0.1),
                                         y_axis_labels = c(0,0.1),
                                         point_size = 4)

THCA_density <- ggdraw(THCA_density) + draw_label(label = "THCA",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/THCA_density.png",plot = THCA_density,base_height = 5/6,base_width = 6.5/4)


# UCEC ---- 

UCEC_density <- density_plotter_function(main_df = combined_all_data.noNA,
                                         this_tumor = "UCEC",
                                         main_size = common.size,
                                         sub_1 = "PTEN R130G",
                                         sub_2 = "NFE2L2 R34G",
                                         sub_1_lab_y = .2,
                                         sub_1_lab_x = 5e5,
                                         sub_1_arrow_x =   259762.83,
                                         sub_1_arrow_y = 0.05,
                                         sub_2_lab_y = .12,
                                         sub_2_lab_x = 5e4,
                                         sub_2_arrow_x =    131243.04,
                                         sub_2_arrow_y = .05,
                                         x_axis_breaks = c(1,1e4,1e5,1e6,1e7),
                                         x_axis_labels = c(expression(10^0),expression(10^4),expression(10^5),expression(10^6),expression(10^7)),
                                         extender_percent = 0.045,
                                         sub_1_arrow_alpha = 1,
                                         sub_2_arrow_alpha = 1,
                                         y_axis_limits = c(0,.4),
                                         y_axis_breaks = c(0,0.1,0.2,0.3),
                                         y_axis_labels = c(0,0.1,0.2,0.3),
                                         point_size = 4)

UCEC_density <- ggdraw(UCEC_density) + draw_label(label = "UCEC",x = .50,y=.9,size = common.size)

save_plot(filename = "figures/UCEC_density.png",plot = UCEC_density,base_height = 5/6,base_width = 6.5/4)


# legend
legend_plot <- ggplot(data = subset(combined_all_data.noNA, tumor_type=="LUAD"), 
                      aes(x=gamma_epistasis)) + 
  geom_density(aes(fill=synonymous),position = "stack") + scale_fill_discrete(name="Substitution\ntype",labels=c("Non-synonymous", "Synonymous")) + theme(legend.title=element_text(size=common.size) , legend.text=element_text(size=common.size))

legend_for_plot <- get_legend(legend_plot)

#### combining density plots  ----

plot_combine <- plot_grid(   UCEC_density,
                             SKCMM_density,
                             COAD_density,
                             STAD_density,
                             SKCMP_density,
                             CESC_density,
                             BRCA_density,
                             BLCA_density,
                             LUAD_density,
                             OV_density,
                           
                             LUSC_density,
                             HNSC_HPVneg_density,
                             LIHC_density,
                             READ_density,
                            
                             GBM_density,
                             
                             LGG_density,
                            
                             
                             ESCA_density,
                             PRAD_density,
                             PAAD_density,
                             
                             KIRC_density,
                             HNSC_HPVpos_density,
                             THCA_density,
                             LAML_density,
                             legend_for_plot,
                             ncol=4,align = 'hv')

combined_with_labs <- ggdraw(plot = plot_grid(NULL,plot_combine,NULL,NULL,nrow = 2,ncol = 2,rel_heights = c(1,0.02),rel_widths = c(0.02,1))) + 
  geom_text(x=.5,y=0.01,label="Selection intensity",size=common.size*(5/14)) + 
  geom_text(x=.01,y=0.5,label="Density",size=common.size*(5/14),angle=90)

save_plot(filename = "figures/combined_density_plot_fill.png",base_width = 6.5,plot = plot_combine,dpi=600,base_height = 10)

save_plot(filename = "figures/combined_density_plot_w_labs_fill.png",base_width = 6.5,plot = combined_with_labs,base_height = 5)




length(which(combined_all_data.noNA$synonymous==T))/nrow(combined_all_data.noNA)

length(which(combined_all_data.noNA$synonymous==F))/nrow(combined_all_data.noNA)
