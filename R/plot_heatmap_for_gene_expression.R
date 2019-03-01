
rm(list=ls())
dat_DEG <- read_delim(pipe("pbpaste"),delim="\t", col_names = TRUE)


# gene     NaCl  `Zn-`    ANM
# AN10628 137.   38.9   87.9
# AN5017   19.0   5.88   8.00
# AN11778  15.5   5.69   6.83
# AN11775  11.0   7.66   7.08
plot_heatmap_for_gene_expression <- function(dat_DEG, output, noClusters, draw_cluster="TRUE"){


          dat_scaled <- dat_DEG %>%
                    gather(Condition,fpkm,-gene) %>%
                    group_by(gene) %>%
                    mutate(zscore=scale(fpkm),Condition=as_factor(Condition)) %>%
                    dplyr::select(-c(fpkm)) %>%
                    spread(Condition, zscore)


          if(draw_cluster=="TRUE"){

                    set.seed(123)

                    dat_scaled <- dat_scaled %>%
                              column_to_rownames("gene") %>%
                              as.matrix()

                    partition = paste0("", kmeans(dat_scaled, centers = noClusters)$cluster)

                    print(table(partition))

                    ht1 <- ComplexHeatmap::Heatmap(dat_scaled,
                                                   name="zscore(RNAP ChIP-seq signal)",
                                                   use_raster = TRUE,
                                                   col = colorRamp2(c( -2,0,2), c("#4d9221","white","#d73027")),
                                                   row_title_gp = gpar(fontsize = 12),
                                                   column_title_rot = 0,
                                                   show_column_names=TRUE,
                                                   show_row_names = FALSE,
                                                   column_names_gp = gpar(fontsize = 12, fontfamily="Arial"),
                                                   heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", fontsize = 12, fontface="bold"),
                                                   cluster_rows = TRUE,
                                                   cluster_columns = FALSE,
                                                   show_row_dend = FALSE)

                    colr <- noClusters+1
                    ht_list = Heatmap(partition, col = structure(2:colr),name="cluster",show_row_names = TRUE, width = unit(3, "mm"))+ht1

                    draw(ht_list,heatmap_legend_side = "bottom", gap = unit(1.5, "mm"), split=partition)

                    pdf(file=paste(output, "withcluster_hm.pdf", sep="_"), width=5, height=9)
                    draw(ht_list,heatmap_legend_side = "bottom", gap = unit(1.5, "mm"), split=partition)
                    dev.off()

                    print("writing clusters ......")
                    dat_DEG$partition = partition
                    write.table(dat_DEG, file=paste(output, "cluster.tab", sep="_"), sep="\t", row.names = FALSE, quote = FALSE)


          }

          else{

                    set.seed(123)
                    km <- dat_scaled %>%
                              column_to_rownames("gene") %>%
                              as.matrix() %>%
                              kmeans(., centers = 4)

                    km_df <- as.tibble(km$cluster) %>%
                              rownames_to_column("geneName")

                    dat_scaled <- dat_scaled %>%
                              left_join(tibble(gene = km_df$geneName , cluster = km_df$value)) %>%
                              arrange(cluster) %>%
                              dplyr::select(-c(cluster)) %>%
                              column_to_rownames("gene") %>%
                              as.matrix()

                    ht1 <- ComplexHeatmap::Heatmap(dat_scaled,
                                                   name="zscore(RNAP ChIP-seq signal)",
                                                   use_raster = TRUE,
                                                   col = colorRamp2(c( -2,0,2), c("#4d9221","white","#d73027")),
                                                   row_title_gp = gpar(fontsize = 12),
                                                   column_title_rot = 0,
                                                   show_column_names=TRUE,
                                                   show_row_names = FALSE,
                                                   column_names_gp = gpar(fontsize = 12, fontfamily="Arial"),
                                                   heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", fontsize = 12, fontface="bold"),
                                                   cluster_rows = TRUE,
                                                   cluster_columns = FALSE,
                                                   show_row_dend = FALSE)

                    draw(ht1,heatmap_legend_side = "bottom", gap = unit(1.5, "mm"))

                    pdf(file=paste(output, "hm.pdf", sep="_"), width=5, height=9)
                    draw(ht1,heatmap_legend_side = "bottom")
                    dev.off()


                    print("writing clusters ......")
                    write.table(km_df, file=paste(output, "cluster.tab", sep="_"), sep="\t", row.names = FALSE, quote = FALSE)

          }
}








