
dat_rna <- tibble::tribble(
             ~Gene,      ~4C,    ~37C,     ~42C,
            "catA", 150.3966, 74.7482,  74.6835,
            "benA",   3.1806,  4.3376,   3.1263,
            "canB",  11.3383, 15.4938,  25.8018,
            "hhtA",  72.6271, 83.8577, 100.0401,
            "hogA",  16.3783, 41.3671,  47.0615,
           "hsp30",   6.5014,  8.4937,  24.2761,
           "hsp70",  71.9213, 79.0689, 119.4281,
           "hsp90",   9.1006, 16.6738,  28.7578,
          "hsp104",  10.4139, 12.0972,  15.2698,
           "hsp20",   7.1911, 34.4908,  61.8152
          )



dat_rna <- read_delim(pipe("pbpaste"),delim="\t", col_names = TRUE)

dat_scale <- dat_rna %>% gather(key=condition, value=fpkm,-Gene) %>%
          mutate(fpkm=ifelse(fpkm< -1,-1,fpkm )) %>%
          # group_by(Gene) %>%
          # mutate(zscore=scale(fpkm),condition=as_factor(condition)) %>%
          # dplyr::select(-c(fpkm)) %>%
          # #arrange(desc(zscore)) %>% # Order high to low
          # ungroup() %>%
          mutate(Gene=as_factor(Gene))

gg=ggplot(data = dat_scale, aes(x = condition, y = Gene)) +
          geom_tile(aes(fill = (fpkm)), colour = "white") +
       #  scale_fill_gradient2(low = "#f46d43",high = "#1a9641") +
          scale_fill_gradient2(low = "#542788",high = "#e08214") +
          theme_bw()+
          scale_y_discrete(position="right")+
          scale_x_discrete(position="top")+
          theme(axis.text.x= element_text(color="black",size=12,family="Arial",angle=0,vjust = 5),
                axis.text.y = element_text(color="black",size=12,family="Arial"),
                axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                legend.position = "bottom",
                legend.title=element_text(color="black",size=12,family="Arial"),
                legend.key.size = unit(1,"line"),
                legend.box.spacing=unit(0.5,"mm"),
                legend.text=element_text(face="bold", color="black",size=10))+
          guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5,
                                        title ="FC(ChIP-seq Signal)", ticks=FALSE))

print(gg)

ggsave(filename = paste("Fig4_A3_heatmap_in_ggplot.pdf"), path="./" ,dpi=300)



ggplot(dat_scale,aes(Gene,fpkm,fill=condition))+geom_col(position="dodge")+
          geom_hline(yintercept = 0)+
          scale_fill_manual(values=c("#377eb8","#d73027"))+
         theme_classic()+
          scale_x_discrete(position="top")+
          theme(axis.text.x= element_text(color="black",size=12,family="Arial",angle=0,vjust=0.8),
                axis.text.y = element_text(color="black",size=12,family="Arial"),
                axis.title.y=element_text(color="black",size=12,family="Arial", face="bold"),
                axis.title.x=element_blank(),
                legend.title=element_blank(),
                legend.key.size = unit(1,"line"),
                legend.position = c(.93, .88),
                legend.box.spacing=unit(0.1,"mm"),
                #legend.position = "bottom",
                legend.text=element_text(color="black",size=12,family="Arial"))+labs(y="log2FC", x="")+
          guides(fill=guide_legend(""))


ggsave(filename = paste("Fig4_A5_heatmap_in_ggplot.pdf"), path="./" ,dpi=300)





