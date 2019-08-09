

dat_scatter <- read_delim(pipe("pbpaste"),delim="\t", col_names = TRUE)

dat_summ <- dat_scatter %>%
          #filter(.[[2]] !=0 & .[[3]]!=0 ) %>%
          gather(sample,fpkm,-geneId) %>%
          separate(col = "sample", into = "replicate",sep = "_spore") %>%
          group_by(geneId, replicate) %>%
          mutate(avrg = (mean(fpkm+0.01))) %>% dplyr::select(-c(fpkm)) %>%
          unique() %>% spread(replicate,avrg)



loadfonts(device = "pdf")

corr_eqn <- function(x,y, digits = 2) {
          corr_coef <- round(cor(x, y), digits = digits)
          paste("italic(r) == ", corr_coef)
}

labels_1 = data.frame(x = 0.1, y =25000, label = corr_eqn((dat_summ$d17_pure), (dat_summ$d17_waterd)))
gg_1 <-   dat_summ %>% ggplot(aes(.[[2]], .[[3]]))+geom_point(size=0.5,alpha=0.5,color= "deepskyblue3")+theme_bw()+
          labs(x="d17_pure spore", y="17d_water spore", title="17d_water")+
          scale_x_log10(breaks=c(1e-01, 1e+01, 1e+03, 1e+05), labels=c("0.1","10","1000","10,000"))+
          scale_y_log10(breaks=c(1e-01, 1e+01, 1e+03, 1e+05), labels=c("0.1","10","1000","10,000"))+
          theme(axis.text.x= element_text(family="Arial", color="black",size=12),
                #axis.ticks = element_blank(),
                axis.text.y = element_text(color="black",size=12,family="Arial"),
                axis.title.y=element_text(color="black",size=12,family="Arial"),
                axis.title.x=element_text(color="black",size=12,family="Arial"),
                legend.title=element_text(color="black",size=12,family="Arial"),
                strip.text = element_text(color="black",size=12,family="Arial"),
                legend.key.size = unit(1.5,"line"),
                plot.title = element_text(hjust = 0.5),
                strip.background = element_rect(colour = "black", fill = "white"),
                legend.text=element_text(color="black",size=12,family="Arial"))+
          guides(fill = guide_legend(title = ""), color=guide_legend(title=""))+
          geom_text(data = labels_1, aes(x = x, y = y,label = label), parse = TRUE,size = 6)
print(gg_1)

ggsave("Fig8_SB_scatterplot.png",dpi=300)

labels_2 = data.frame(x = 0.1, y =25000, label = corr_eqn((dat_summ$d3), (dat_summ$d17_pure)))
gg_2 <-   ggplot(dat_summ,aes(d3, d17_pure))+geom_point(alpha=0.5,color= "deepskyblue3")+theme_bw()+
          labs(x="3d spore", y="17d_dispersal spore", title="17d_dispersal")+
          scale_x_log10(breaks=c(1e-01, 1e+01, 1e+03, 1e+05), labels=c("0.1","10","1000","10,000"))+
          scale_y_log10(breaks=c(1e-01, 1e+01, 1e+03, 1e+05), labels=c("0.1","10","1000","10,000"))+
          theme(axis.text.x= element_text(family="Arial", color="black",size=10),
                #axis.ticks = element_blank(),
                axis.text.y = element_text(color="black",size=12,family="Arial"),
                axis.title.y=element_text(color="black",size=12,family="Arial"),
                axis.title.x=element_text(color="black",size=12,family="Arial"),
                legend.title=element_text(color="black",size=12,family="Arial"),
                strip.text = element_text(color="black",size=12,family="Arial"),
                legend.key.size = unit(1.5,"line"),
                plot.title = element_text(hjust = 0.5),
                strip.background = element_rect(colour = "black", fill = "white"),
                legend.text=element_text(color="black",size=12,family="Arial"))+
          guides(fill = guide_legend(title = ""), color=guide_legend(title=""))+
          geom_text(data = labels_2, aes(x = x, y = y,label = label), parse = TRUE,size = 6)
print(gg_2)

ggsave("Fig8_B_scatterplot.png",dpi=300)



dd <- dat_scatter %>% column_to_rownames(geneID)


GGally::ggpairs(dat_scatter[,2:7])





