
#---- Line plot
gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"

mylist <- read_delim(pipe("pbpaste"), delim="\t", col_names=FALSE)

#----


lineplot_for_bw <- function(gff_file, mylist,bw_file, output, tss="TRUE"){
          library(tidyverse)

          gff <- makeTxDbFromGFF(gff_file, metadata = T)
          genes <- genes(gff)

          #--- provide the genelist data with expression value for the list one

          mylist <- mylist %>% distinct(.) %>% arrange(desc(X2))
          genes_1 <-  subset(genes, genes$gene_id %in% mylist$X1)
          genes_1 <- genes_1[match(mylist$X1,genes_1$gene_id),]

          mylist <- mylist %>% mutate(quantile=ntile(X2,4))

          genes_1$quantile = mylist$quantile[match(genes_1$gene_id,mylist$X1)]

          if(tss=="TRUE"){
                    tss <- promoters(genes_1, upstream = 0, downstream = 1)
                    breaks=c("u2","u50","d50")
                    labels=c("-1kb","TSS","+1kb")
                    inter=50
          }
          else{
                    tss=genes_1
                    breaks=c("u2","u50","d2","d50")
                    labels=c("-1kb","TSS","TES","+1kb")
                    inter =c(50,117)
          }

          bw <- import(bw_file)

          mat_1 <- EnrichedHeatmap::normalizeToMatrix(bw, tss, value_column = "score",background = 0,
                                                      smooth = TRUE,extend = c(1000))

          print(mat_1)

          norm_mat <- data.frame(mat_1)

          summarised_mat <- norm_mat %>%
                    rownames_to_column("geneName") %>% as_tibble() %>%
                    left_join(tibble(geneName = genes_1$gene_id , quantile = genes_1$quantile)) %>%
                    #mutate(quantile = genes_1$quantile[match(row.names(norm_mat), genes_1$gene_id)]) %>% as_tibble() %>%
                    gather(key=column_name, value=value,-quantile,-geneName) %>% mutate(column_name=as_factor(column_name)) %>%
                    group_by(column_name,quantile) %>%
                    summarise(median=median(value)) %>% ungroup()



          print(summarised_mat)

          gg <- ggplot(summarised_mat,aes(column_name,median, group=quantile,color =rev(quantile)))+
                    geom_line(lwd=1.1)+
                   # scale_color_gradient(high="#fd8d3c", low="#b10026")+
                    #scale_color_gradient(low="#7a0177", high="#fcc5c0")+
                    #scale_color_gradient(high="#b2e2e2", low="#005824")+
                    #scale_color_gradient(high="#e5f5e0", low="#005a32")+
                    scale_x_discrete(breaks=breaks,labels=labels,"")+
                    theme_bw()+
                    geom_vline(xintercept = inter, size=0.8,color="red",linetype = "dashed")+
                    theme(
                              panel.grid.minor.y = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.x= element_blank(),
                              axis.text.x= element_text(color="black",size=12,family="Arial",angle=0,vjust=0.8),
                              axis.text.y = element_text(color="black",size=12,family="Arial"),
                              axis.title.y=element_text(color="black",size=12,family="Arial", face="bold"),
                              axis.title.x=element_blank(),
                              legend.title=element_text(color="black",size=12,family="Arial"),
                              legend.key.size = unit(1,"line"),
                              legend.box.spacing=unit(0.5,"mm"),
                              legend.text=element_text(color="black",size=12,family="Arial")) +
                    ylab("median (TBP binding)")+
                    guides(color=guide_legend(title=""))

          print(gg)

          ggsave(filename = paste(output,"_lineplot.pdf",sep=""), device = "pdf", path="./", width = 10, height =8,unit="cm", dpi=300)

}

