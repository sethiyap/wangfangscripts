####
# This function allows you to plot multiple genelists average from same bw file
# Please specify the gff file and load your genelist with groups
#---
gff_file <- "Org_version_s10-m04-r07_features.gff"
mylist <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE)
#Input list format
# AN5056  group_1
# AN4954  group_1
# AN4871  group_1
# AN4848  group_2
# AN4770  group_2
# AN4716  group_3

#Provide bw_file path
bw_file <- "sample_normalized.bw"
#-----
#Run function as
multiple_genelist_line_from_same_bw(gff_file = gff_file,mylist = mylist, bw_file = bw_file ,output = "Test", tss = FALSE,mean = TRUE)

#Load function first
multiple_genelist_line_from_same_bw <- function(gff_file, mylist,bw_file, output, tss="TRUE", mean="TRUE"){

          library(tidyverse)
          library(EnrichedHeatmap)
          library(GenomicFeatures)
          library(rtracklayer)

          gff <- makeTxDbFromGFF(gff_file, metadata = T, format="gff3")
          genes <- genes(gff)

          #--- provide the genelist data with expression value for the list one

          mylist <- mylist %>% distinct(.)
          genes_1 <-  subset(genes, genes$gene_id %in% mylist$X1)

          genes_1$group = mylist$X2[match(genes_1$gene_id,mylist$X1)]

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

          if(mean=="TRUE"){

          print("Computing mean...")

          summarised_mat <- norm_mat %>%
                    rownames_to_column("geneName") %>% as_tibble() %>%
                    left_join(tibble(geneName = genes_1$gene_id , group = genes_1$group)) %>%
                    gather(key=column_name, value=value,-group,-geneName) %>% mutate(column_name=as_factor(column_name)) %>%
                    group_by(column_name,group) %>%
                    summarise(average=mean(value)) %>% ungroup()

          print(summarised_mat)

          gg <- ggplot(summarised_mat,aes(column_name,average, group=group,color =rev(group)))+
                    geom_line(lwd=1.1)+
                    scale_x_discrete(breaks=breaks,labels=labels,"")+
                    theme_bw()+
                    geom_vline(xintercept = inter, size=0.8,color="red",linetype = "dashed")+
                    theme(
                              panel.grid.minor.y = element_blank(),
                              panel.grid.major.y = element_blank(),
                              panel.grid.minor.x = element_blank(),
                              panel.grid.major.x= element_blank(),
                              axis.text.x= element_text(color="black",size=12,angle=0,vjust=0.8),
                              axis.text.y = element_text(color="black",size=12),
                              axis.title.y=element_text(color="black",size=12, face="bold"),
                              axis.title.x=element_blank(),
                              legend.title=element_text(color="black",size=12),
                              legend.key.size = unit(1,"line"),
                              legend.box.spacing=unit(0.5,"mm"),
                              legend.text=element_text(color="black",size=12)) +
                    ylab("mean signal")+
                    guides(color=guide_legend(title=""))
          }

          else{
                    print("Computing median...")
                    summarised_mat <- norm_mat %>%
                              rownames_to_column("geneName") %>% as_tibble() %>%
                              left_join(tibble(geneName = genes_1$gene_id , group = genes_1$group)) %>%
                              gather(key=column_name, value=value,-group,-geneName) %>% mutate(column_name=as_factor(column_name)) %>%
                              group_by(column_name,group) %>%
                              summarise(average=median(value)) %>% ungroup()



                    print(summarised_mat)

                    gg <- ggplot(summarised_mat,aes(column_name,average, group=group,color =rev(group)))+
                              geom_line(lwd=1.1)+
                              scale_x_discrete(breaks=breaks,labels=labels,"")+
                              theme_bw()+
                              geom_vline(xintercept = inter, size=0.8,color="red",linetype = "dashed")+
                              theme(
                                        panel.grid.minor.y = element_blank(),
                                        panel.grid.major.y = element_blank(),
                                        panel.grid.minor.x = element_blank(),
                                        panel.grid.major.x= element_blank(),
                                        axis.text.x= element_text(color="black",size=12,angle=0,vjust=0.8),
                                        axis.text.y = element_text(color="black",size=12),
                                        axis.title.y=element_text(color="black",size=12, face="bold"),
                                        axis.title.x=element_blank(),
                                        legend.title=element_text(color="black",size=12),
                                        legend.key.size = unit(1,"line"),
                                        legend.box.spacing=unit(0.5,"mm"),
                                        legend.text=element_text(color="black",size=12)) +
                              ylab("median signal")+
                              guides(color=guide_legend(title=""))
          }
          print(gg)

          ggsave(filename = paste(output,"_lineplot.pdf",sep=""), device = "pdf", path="./", width = 10, height =8,unit="cm", dpi=300)


}
