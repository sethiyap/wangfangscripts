#-- bubble plot ----
data <- read_delim(pipe("pbpaste"), delim="\t", col_names =T)

plot_go_in_diamond <- function(data,output_name){

          #--- Load libraries
          library(ggplot2)
          library(tidyverse)
          library(extrafont)
          loadfonts(device = "pdf")

          print(head(data))


          data_melt <- data %>%
                    mutate(Percent.Enrichment=100*(data$Significant/data$Annotated))

if("Class" %in% colnames(data)){
          data_melt$Term <- factor(data$Term, levels=rev(unique(data$Term)))
          data_melt$Category <- factor(data$Category, levels=(unique(data$Category)))
          data_melt$Class<- factor(data$Class, levels=(unique(data$Class)))

          gg= ggplot(data_melt , aes(y=Term, x=Class,ordered=TRUE, fill=Category))+xlab("")+ylab("")+
                    geom_point(aes(size=Percent.Enrichment,alpha=-log10(Pvalue)),shape=23)+
                    scale_size(range = c(1,15))+
                    #scale_fill_manual(values=c("#c51b7d","#01665e"))+
                    scale_fill_manual(values=c("#e41a1c","#1b9e77","#f0027f"))+
                    scale_y_discrete(position="left")+
                    theme_bw()+
                    scale_x_discrete(position="top")+
                    theme(axis.text.x= element_text(color="black",size=12,family="Arial"),
                          axis.text.y = element_text(color="black",size=12,family="Arial"),
                          axis.title.x=element_text(color="black",size=12,family="Arial"),
                          legend.title=element_text(color="black",size=12,family="Arial"),
                          legend.key.size = unit(0.5,"line"),
                          panel.spacing = unit(1, "lines"),
                          legend.text=element_text(color="black",size=12,family="Arial"))+
                    scale_color_manual(values ="black")+
                    guides(fill = guide_legend(title="",override.aes = list(size=8)),
                           size=guide_legend(title="percentage of \ngene over bgd",override.aes = list(size=c(2,4,6,8))),
                           alpha=guide_legend(override.aes=list(size=5)))
          print(gg)
          ggsave(paste(output_name,"GoByColorPvalue.pdf",sep=""), height = 5.5, width=8)

          }
          else{
                    data_melt$Term <- factor(data$Term, levels=rev(unique(data$Term)))
                    data_melt$Category <- factor(data$Category, levels=(unique(data$Category)))

                    gg= ggplot(data_melt , aes(y=Term, x=Category,ordered=TRUE, fill=Category))+xlab("")+ylab("")+
                              geom_point(aes(size=Percent.Enrichment,alpha=-log10(Pvalue)),shape=23)+
                              scale_size(range = c(1,15))+
                              scale_fill_manual(values=c("#e41a1c","#1b9e77","#f0027f"))+
                              scale_y_discrete(position="right")+
                              theme_bw()+
                              scale_x_discrete(position="top")+
                              theme(axis.text.x= element_text(color="black",size=12,family="Arial"),
                                    axis.text.y = element_text(color="black",size=12,family="Arial"),
                                    axis.title.x=element_text(color="black",size=12,family="Arial"),
                                    legend.title=element_text(color="black",size=12,family="Arial"),
                                    legend.key.size = unit(0.5,"line"),
                                    panel.spacing = unit(1, "lines"),
                                    legend.text=element_text(color="black",size=12,family="Arial"))+
                              scale_color_manual(values ="black")+
                              guides(fill = guide_legend(title="",override.aes = list(size=8)),
                                     size=guide_legend(title="percentage of \ngene over bgd",override.aes = list(size=c(2,4,6,8))),
                                     alpha=guide_legend(override.aes=list(size=5)))
                    print(gg)
                    ggsave(paste(output_name,"GoByColorPvalue.pdf",sep=""), height = 10, width=9)

          }
}


tibble::tribble(
                                                               ~Term, ~Annotated, ~Significant, ~Expected,  ~Pvalue, ~Category, ~Class,
                                                       "translation",        457,           55,     10.52,  4.2e-26,      "Zn", "DOWN",
          "maturation of SSU-rRNA from tricistronic rRNA transcript",         52,           17,       1.2,  2.9e-14,      "Zn", "DOWN",
                                        "ribosomal subunit assembly",         11,            6,      0.25,  5.9e-08,      "Zn", "DOWN",
                                          "rRNA export from nucleus",         11,            6,      0.25,  5.9e-08,      "Zn", "DOWN",
                    "formation of translation preinitiation complex",          5,            2,      0.12,  0.00504,      "Zn", "DOWN",
                                                  "ncRNA processing",        205,           25,      4.72,  0.01774,      "Zn", "DOWN",
                                                     "ADP transport",          1,            1,      0.02,  0.02302,      "Zn", "DOWN",
                                                   "tRNA processing",         76,            5,      1.75,  0.03059,      "Zn", "DOWN",
                                        "response to osmotic stress",         80,            5,      1.84,    0.037,      "Zn", "DOWN",
                    "positive regulation of protein kinase activity",         15,            2,      0.35,  0.04551,      "Zn", "DOWN",
                            "zinc ion import across plasma membrane",          2,            2,      0.02,  6.4e-05,      "Zn",   "UP",
                                                     "intron homing",          3,            2,      0.02,  0.00019,      "Zn",   "UP",
                                              "response to zinc ion",          5,            2,      0.04,  0.00063,      "Zn",   "UP",
                                           "amino acid biosynthesis",        138,            4,      1.11,   0.0253,      "Zn",   "UP"
          )
