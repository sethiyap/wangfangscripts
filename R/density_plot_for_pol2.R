Pol2_dat <- read_delim(pipe("pbpaste"),delim="\t", col_names = TRUE)

#' Title density_plot_for_pol2
#'
#' @param Pol2_dat add Pol2 data in tabular format with Gene, sample_1, sample_2
#' @param outfile pdf of the desnity plot
#' @param background_value determine the background to highlight the background value
#'
#' @return
#' @export
#' @import tidyverse
#' @import extrafont
#' @author pooja sethiya
#' @examples
#'
density_plot_for_pol2 <- function(Pol2_dat, outfile, background_value, title){

          library(tidyverse)
          library(extrafont)
          loadfonts(device = "pdf")

          #--- rename the coumns
          colnames(Pol2_dat) <- c("Gene","RNAP Signal", "Input")

          #--- reshape data for plotting
          ppm <- Pol2_dat %>% gather(key=variable, value=value,-Gene)

          #--- plot ggplot
          gg = ggplot(ppm,aes(value+0.01,fill=variable))+geom_density(alpha=0.5)+
                    xlab("RNAP ChIP-seq signal")+
                    ylab("Density of genes")+
                    theme_bw()+
                    scale_fill_manual(values=c("#ffeda0","#f03b20"))+
                    ggtitle(label = title)+
                    geom_vline(xintercept=background_value,size = 1, colour = c("grey"),linetype = "dashed",show.legend = TRUE)+
                    scale_x_sqrt(expand = c(0,0), limits = c(0,max(ppm$value)))+
                    theme(axis.text.x= element_text( colour="black", size=12,angle=0,vjust=0.8,family="Arial"),
                          axis.text.y = element_text( color="black",size=12,family="Arial"),
                          axis.title.y=element_text(color="black",size=12,family="Arial"),
                          axis.title.x=element_text(color="black",size=12,family="Arial"),
                          plot.title = element_text(hjust = 0.5,color="black",size=12,family="Arial", face="italic"),
                          legend.title=element_blank(),
                          legend.key.size = unit(1,"line"),
                          legend.position = c(.80, .85),
                          legend.text=element_text( color="black",size=12,family="Arial"))+
                    geom_text(aes(x=background_value,label=paste("\n background =", round(background_value,2))), y=0.4,angle=90,colour="black", size=6,family="Arial")+
                    guides(fill = guide_legend(title=""))
          print(gg)

          #--- save file
          ggsave(filename = paste(outfile,"_Pol2_distribution.pdf",sep=""), path="./",width=5, height=5, dpi=300)
}



# density_plot_for_pol2(Pol2_dat,"Tm_Pol2",9.1,"Tm")





