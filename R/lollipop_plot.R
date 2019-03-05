


dat <- read_delim(pipe("pbpaste"),delim="\t", col_names = TRUE)

lollipop_plot <- function(dat, output_name, title){

          library(tidyverse)

          #--- process data
          dat_plot <-   dat %>% gather(key="Sample", value="fpkm",-Gene) %>%
                    group_by(Gene) %>%
                    mutate( mymean = mean(fpkm)) %>% ungroup() %>%
                    arrange((mymean)) %>%  # Order high to low
                    mutate(Gene=factor(Gene, unique(Gene)))


          # plot
         gg <-  ggplot(dat_plot,aes(x = Gene, y = fpkm, color=Sample)) +
                    geom_line(aes(group=Gene), color="skyblue", size=0.8) +
                    geom_point(size=4, alpha=0.8) +
                    scale_color_manual(values=c("#1b7837","#762a83"))+
                    labs(y="RNAP ChIP-seq Signal", title=title)+
                    coord_flip() +
                    theme_bw() +
                    theme(panel.border = element_blank(),
                          axis.text.x= element_text(color="black",size=12,family="Arial",angle=0,vjust=0.8),
                          axis.text.y = element_text(color="black",size=12,family="Arial"),
                          axis.title.y=element_blank(),
                          plot.title = element_text(hjust = 0.5,color="black",size=12,family="Arial"),
                          legend.title=element_blank(),
                          legend.key.size = unit(1,"line"),
                          legend.box.spacing=unit(0.5,"mm"),
                          legend.position = c(0.89, 0.1),
                          legend.text=element_text(color="black",size=12,family="Arial"))+
                    guides(color = guide_legend(title=""))

         print(gg)

         ggsave(filename = paste(output_name,"_lollipop_in_ggplot.pdf",sep=""), path="./" ,dpi=300)


}



