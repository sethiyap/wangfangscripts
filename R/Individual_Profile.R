
gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"

mylist <- read_delim(pipe("pbpaste"), delim="\t", col_names=FALSE)

bw_files_dir <- (".")

bw_files <- list.files(bw_files_dir, pattern = "Spore_.*_normalized.*", recursive = T, full.names = T)

#--- replace file names to get the sample names
#eg: ./H2A.Z_by_Htz1_NH4_GTTGTCCCA_CW466_525_normalized.bw
# print only H2A.Z_by_Htz1_NH4

names(bw_files) <- gsub(pattern = "_[[:upper:]]{6,}_.*_normalized.bw",replacement = "", bw_files)
names(bw_files) <- gsub(pattern = "./",replacement = "", names(bw_files))

#names(bw_files) <- c("hepA_HA_20hr","hepA_HA_48hr","hypha_h3k9me3_wt","hypha_h3k9me3_hepA_del", "spore_h3k9me3_wt","spore_h3k9me3_hepA_del" )
bw_files <- broom::tidy(bw_files)
print(bw_files)

gff <- makeTxDbFromGFF(gff_file, metadata = T)
genes <- genes(gff)

#--- provide the genelist data with expression value for the list one

mylist <- mylist %>% arrange(desc(X2))
genes_1 <-  subset(genes, genes$gene_id %in% mylist$X1)
genes_1 <- genes_1[match(mylist$X1,genes_1$gene_id),]

#----- normalise mat
tss = promoters(genes_1, upstream = 0, downstream = 1)

xx <- bw_files %>%
          dplyr::mutate(bw = map(x, function(ii) {
                    rtracklayer::import(ii)
          })) %>%
          dplyr::mutate(norm_matrix = purrr::map(bw, function(ii) {
                    nn <- EnrichedHeatmap::normalizeToMatrix(ii, tss, value_column = "score",background = 0,
                                                             smooth = TRUE,extend = c(1000))
                    nn[nn<0]=0
                    return(nn)
          }))



eml_1 <- EnrichedHeatmap::EnrichedHeatmap(log2(xx$norm_matrix[[1]]+0.01),
                                name = xx$names[[1]],
                                column_title = xx$names[[1]],
                                row_order = order(mylist$X2, decreasing = TRUE),
                         cluster_rows = FALSE,
                         show_row_names = FALSE,
                         axis_name_rot = 90,
                         heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", legend_width = unit(3, "cm"),
                                                     title_position = "topcenter",labels_gp = gpar(fonsize=12, fontfamily="Arial")),
                         axis_name = c("-1kb","TSS", "+1kb"),
                         axis_name_gp = gpar(fonsize=12, fontfamily="Arial"),
                         col = colorRamp2(breaks = c(0,4,6,8,10),
                                          colors = c("white","#fef0d9","#ef6548","#d7301f","#990000")),
                         top_annotation = HeatmapAnnotation(lines = anno_enriched(axis_param =list( facing="inside",side="left",gp=gpar(fonsize=12, fontfamily="Arial")),
                                                                                  ylim = c(2.8,5.5),height = unit(2, "cm")
                                                                                 ))
                         )





eml_2 <-  EnrichedHeatmap::EnrichedHeatmap(log2(xx$norm_matrix[[2]]+0.01),
                                           name = xx$names[[2]],
                                           column_title = xx$names[[2]],
                                                   cluster_rows = FALSE,
                                                   show_row_names = FALSE,
                                                   axis_name_rot = 90,
                                                   heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", legend_width = unit(3, "cm"),
                                                                               title_position = "topcenter",labels_gp = gpar(fonsize=12, fontfamily="Arial")),
                                                   axis_name = c("-1kb","TSS", "+1kb"),
                                                   axis_name_gp = gpar(fonsize=12, fontfamily="Arial"),
                                           col = colorRamp2(breaks = c(0,1,2,3,4),
                                                                 colors = c("white","#b2e2e2","#99d8c9","#66c2a4", "#005824")),
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(axis_param =list( facing="inside",side="left",gp=gpar(fonsize=12, fontfamily="Arial")),
                                                                                           ylim = c(2.8,5.5),height = unit(2, "cm")
                                                                                                  ))
                                                  )


eml_3 <- EnrichedHeatmap::EnrichedHeatmap(log2(xx$norm_matrix[[3]]+0.01),
                                          name = xx$names[[3]],
                                          column_title = xx$names[[3]],
                                                   cluster_rows = FALSE,
                                                   show_row_names = FALSE,
                                                   axis_name_rot = 90,
                                                   heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", legend_width = unit(3, "cm"),
                                                                               title_position = "topcenter",labels_gp = gpar(fonsize=12, fontfamily="Arial")),
                                                   axis_name = c("-1kb","TSS", "+1kb"),
                                                   axis_name_gp = gpar(fonsize=12, fontfamily="Arial"),
                                         col = colorRamp2(breaks = c(1,2,3,4,5),
                                                          colors = c("white","#d0d1e6","#a6bddb","#0570b0","#034e7b")),
                                         top_annotation = HeatmapAnnotation(lines = anno_enriched(axis_param =list( facing="inside",side="left",gp=gpar(fonsize=12, fontfamily="Arial")),
                                                                                                  ylim = c(2.8,5.5),height = unit(2, "cm")
                                                                                             ))
                                         )

eml_list <- eml_1+eml_2+eml_3


print("plotting....")
pdf(file=paste("FIG2_B3", length(genes_1), "hm.pdf", sep="_"), width=3, height=10)
draw(eml_3, heatmap_legend_side = "top", gap = unit(2, "mm"))
dev.off()

