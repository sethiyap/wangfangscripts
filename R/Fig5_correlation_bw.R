library(reshape2)
#--- get correlation matrix by computing correlation between bw files using bigWigCorrelate command on shell

dat_correl <- read_delim(pipe("pbpaste"),delim="\t", col_names=TRUE)


# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
          cormat[upper.tri(cormat)] <- NA
          return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
          cormat[lower.tri(cormat)]<- NA
          return(cormat)
}

upper_tri <- dat_correl %>% column_to_rownames("Time") %>% as.matrix() %>% get_upper_tri(.)

melted_cormat <- melt(upper_tri, na.rm = TRUE)



gg <- ggplot(melted_cormat,aes(Var2,Var1, fill=value))+
          geom_tile(color = "white")+
          scale_fill_gradient2(low = "#4393c3", high = "#d73027", limit = c(0,1), space = "Lab",
                               name="Pearson\nCorrelation") +
          theme_minimal()+ # minimal theme
          theme(axis.text.x = element_text(angle = 45, vjust=1,hjust=1,size = 12,color="black",family="Arial", face="bold"),
                axis.text.y = element_text(angle = 0, vjust=1,hjust=1,size = 12,color="black",family="Arial", face="bold"))+
          coord_fixed()



gg +
          geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 4,family="Arial") +
          theme(    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),
                    axis.ticks = element_blank(),
                    legend.justification = c(1, 0),
                    legend.position = c(0.4, 0.7),
                    legend.direction = "horizontal")+
          guides(fill = guide_colorbar(barwidth = 7, barheight = 1,ticks=FALSE,
                                       title.position = "top", title.hjust = 0.5))


ggsave("Fig5_genes_correlation.pdf",dpi=300, width=5, height = 5)



