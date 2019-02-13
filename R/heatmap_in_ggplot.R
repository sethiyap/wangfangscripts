heatmap_in_ggplot <-
function(dat, output_name){

          #--- Load libraries
          library(ggplot2)
          library(tidyverse)

          #--- Data format
          # Gene	Spore	Hypha	Input	Category
          # cetL	111.7443	2.8244	5.3759	Conidia enriched transcript
          # wetA	42.2475	3.1493	6.7317	Spore maturation
          # wA	36.6445	2.2049	4.3408	Spore maturation

          #--- facet if different Categories provided
          if ("Category" %in% colnames(dat)){

                    message("plotting categories...")
                    #--- Gather data for plotting
                    dat_plot <-   dat %>% gather(key="Sample", value="fpkm",-Gene, -Category) %>%
                              mutate(Sample=as_factor(Sample), Category=as_factor(Category)) %>% group_by(Category,Sample) %>%
                              arrange((fpkm),.by_group = TRUE)

                    dat_plot$Gene <- factor(dat_plot$Gene, levels=unique(dat_plot$Gene))
                    #--- Plot heatmap
                   gg= ggplot(data = dat_plot, aes(x = Sample, y = Gene)) +
                              geom_tile(aes(fill = log2(fpkm)), colour = "white") +
                              facet_wrap(~Category, nrow=length(unique(dat_plot$Category)),scales = "free_y", strip.position="left") +
                              scale_fill_gradient(low = "#ffeda0",high = "#f03b20", limits = c(0,max(log2(dat_plot$fpkm)))) +theme_bw()+
                              scale_y_discrete(position="right")+
                              theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=0,vjust=0.8),
                                    axis.text.y = element_text(face="bold", color="black",size=10),
                                    axis.title.y=element_blank(),
                                    axis.title.x=element_blank(),
                                    legend.position = "top",
                                    #strip.background = element_rect( fill="white"),
                                    strip.text=element_text(face="bold", color="black",size=11),
                                    legend.title=element_text(face="bold", color="black",size=10),
                                    legend.key.size = unit(1,"line"),
                                    legend.box.spacing=unit(0.5,"mm"),
                                    legend.text=element_text(face="bold", color="black",size=10))+
                              guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,title ="log2(RNA Pol-II Signal)"))

                   print(gg)
                   ggsave(filename = paste(output_name,"_heatmap_in_ggplot.png"), device = "png", path="./", width = 8, height =19,unit="cm", dpi=300)

          }

          else{
                    message("no category provided! plotting data...")
                    dat_plot <-   dat %>% gather(key="Sample", value="fpkm",-Gene) %>%
                                        mutate(Sample=as_factor(Sample)) %>% group_by(Sample) %>% arrange((fpkm),.by_group=TRUE) # Order high to low

                    dat_plot$Gene <- factor(dat_plot$Gene, levels=unique(dat_plot$Gene))
                  gg=ggplot(data = dat_plot, aes(x = Sample, y = Gene)) +
                              geom_tile(aes(fill = log2(fpkm)), colour = "white") +
                              scale_fill_gradient(low = "#ffeda0",high = "#f03b20",limits = c(0,max(log2(dat_plot$fpkm)))) +theme_bw()+
                              scale_y_discrete(position="right")+
                              theme(axis.text.x= element_text(face="bold", colour="black", size=10,angle=0,vjust=0.8),
                                    axis.text.y = element_text(face="bold", color="black",size=10),
                                    axis.title.y=element_blank(),
                                    axis.title.x=element_blank(),
                                    legend.position = "top",
                                    legend.title=element_text(face="bold", color="black",size=10),
                                    legend.key.size = unit(1,"line"),
                                    legend.box.spacing=unit(0.5,"mm"),
                                    legend.text=element_text(face="bold", color="black",size=10))+
                              guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,title ="log2(RNA Pol-II Signal)"))

                     print(gg)

                     ggsave(filename = paste(output_name,"_heatmap_in_ggplot.png"), device = "png", path="./", width = 8, height =8,unit="cm", dpi=300)

          }



}
