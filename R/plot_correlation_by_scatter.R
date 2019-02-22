library(tidyverse)

#supp correlation ----
dat_scatter <- read_delim(pipe("pbpaste"), delim="\t", col_names = TRUE)

# Gene     Set1  Set2
# AN7895   3.43  4.88
# AN8366   2.43  4.60
# AN0833   8.74  9.27
# AN11702  3.07  6.40

#' plot_correlation_by_scatter
#'
#' @param dat_scatter provide data in tidy format, and mention category for multi sample scatterplot
#' @param output_name provide output file name
#' @param plot_title title for the plot
#' @param color color for the points to be plotted
#' @import tidyverse
#' @return
#' @author pooja sethiya
#' @export
#'
#' @examples
plot_correlation_by_scatter <- function(dat_scatter, output_name, plot_title="Sample", color="black"){

          colname=colnames(data)
          corr_eqn <- function(x,y, digits = 2) {
                    corr_coef <- round(cor(x, y), digits = digits)
                    paste("italic(r) == ", corr_coef)
          }

          if ("Category" %in% colnames(dat_scatter)){
                    gg <- ggplot(dat_scatter,aes(log2(Set1), log2(Set2)))+geom_point(alpha=0.5, color=color)+
                              facet_wrap(~Class)+
                              labs(x="Set1", y="Set2", title=plot_title)+
                              theme_bw()+
                              theme(axis.text.x= element_text(face="bold", color="black",size=10),
                                    axis.text.y = element_text(face="bold", color="black",size=10),
                                    axis.title.y=element_text(face="bold", color="black",size=10),
                                    axis.title.x=element_text(face="bold", color="black",size=10),
                                    legend.title=element_text(face="bold", color="black",size=10),
                                    strip.text = element_text(face="bold", color="black",size=10),
                                    legend.key.size = unit(1,"line"),
                                    plot.title = element_text(hjust = 0.5),
                                    strip.background = element_rect(colour = "black", fill = "white"),
                                    legend.text=element_text(face="bold", color="black",size=14))+
                              guides(fill = guide_legend(title = ""), color=guide_legend(title=""))
                    print(gg)
                    ggsave(filename = paste(output_name,"_correlationByScatter.pdf"), path="./",unit="cm", dpi=300)

          }

          else{
                    labels = data.frame(x = 1, y = 8, label = corr_eqn(log2(dat_scatter$Set1), log2(dat_scatter$Set2)))
                    gg <- ggplot(dat_scatter,aes(log2(Set1), log2(Set2)))+geom_point(alpha=0.5, color=color)+
                              theme_bw()+labs(x="Set1", y="Set2", title=plot_title)+
                              theme(axis.text.x= element_text(face="bold", color="black",size=10),
                                    #axis.ticks = element_blank(),
                                    axis.text.y = element_text(face="bold", color="black",size=10),
                                    axis.title.y=element_text(face="bold", color="black",size=10),
                                    axis.title.x=element_text(face="bold", color="black",size=10),
                                    legend.title=element_text(face="bold", color="black",size=10),
                                    strip.text = element_text(face="bold", color="black",size=10),
                                    legend.key.size = unit(1.5,"line"),
                                    plot.title = element_text(hjust = 0.5),
                                    strip.background = element_rect(colour = "black", fill = "white"),
                                    legend.text=element_text(face="bold", color="black",size=14))+
                              guides(fill = guide_legend(title = ""), color=guide_legend(title=""))+
                              geom_text(data = labels, aes(x = x, y = y,label = label), parse = TRUE,size = 6)
                    print(gg)

                    ggsave(filename = paste(output_name,"_correlationByScatter.pdf"), path="./",unit="cm", dpi=300)

          }
}
