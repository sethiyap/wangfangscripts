gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"
mylist <- read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE)
bw_files_dir <-c(".")


#--- packages
library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(rbamtools)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)
library(extrafont)
loadfonts(device = "pdf")

gff <- makeTxDbFromGFF(gff_file, metadata = T)
genes <- genes(gff)

          if(mylist== "NULL"){
                    genes <- genes
          }

          else {
                    genes <-  subset(genes, genes$gene_id %in% mylist$X1)

                    genes <- genes[match(mylist$X1,genes$gene_id),]

                    row_order = order(mylist$X2, decreasing=TRUE)
          }


tss = promoters(genes, upstream = 0, downstream = 1)
## prepare signal data
bw_files <- list.files(bw_files_dir, pattern = "*_normalized.*", recursive = T, full.names = T)

#--- replace file names to get the sample names
#eg: ./H2A.Z_by_Htz1_NH4_GTTGTCCCA_CW466_525_normalized.bw
# print only H2A.Z_by_Htz1_NH4

names(bw_files) <- c("d3_RNA", "d17_RNA","3d_Pol2","17d_Pol2","30d_Pol2")

bw_files <- broom::tidy(bw_files)
print(bw_files)



## generate normalised matrix in tidy way
xx <- bw_files %>%
          dplyr::mutate(bw = map(x, function(ii) {
                    rtracklayer::import(ii)
          })) %>%
          dplyr::mutate(norm_matrix = purrr::map(bw, function(ii) {
                    nn <- EnrichedHeatmap::normalizeToMatrix(ii, tss, value_column = "score", mean_mode = "w0",
                                                             keep = c(0.0, 0.99),background = 0,
                                                             smooth = TRUE,extend = c(1000, 3000), w = 10)
                    nn[nn<0]=0
                    return(nn)
          }))


xx_1 <- xx[c(1:2),]

xx_2 <- xx[c(3:5),]
get_enrichment_heatmap_list <- function(x, names, titles, ...) {
          #x <- xx$norm_matrix

          ll <- length(x)

          ## first heatmap
          ehml <- EnrichedHeatmap(mat = x[[1]], name = names[[1]], column_title = titles[[1]], show_heatmap_legend = T,
                                  use_raster = TRUE,
                                  # col=colorRamp2(quantile(x[[1]], c(0.00,0.25,0.50,0.75,0.9,0.99)),
                                  #                col=c("white","#dadaeb","#bcbddc", "#9e9ac8","#807dba", "#3f007d")),
                                  col=colorRamp2(quantile(x[[1]], c(0.00,0.25,0.50,0.75,0.9,0.99)),
                                                 col=c("#c994c7","#fde0dd","#fcc5c0", "#fa9fb5", "#ae017e","#7a0177")), #for pol2
                                   ...)

          ## several other heatmaps if length of x > 1.
          if (ll > 1) {
                    for (i in 2:ll) {
                              print(i)
                              ehml <- ehml +
                                        EnrichedHeatmap(
                                                  mat = x[[i]],
                                                  name = ifelse(length(names) >= i, names[i], "NA"),
                                                  use_raster = TRUE,
                                                  # col=colorRamp2(quantile(x[[i]], c(0.00,0.25,0.50,0.75,0.9,0.99)),
                                                  #                col=c("white","#dadaeb","#bcbddc", "#9e9ac8","#807dba", "#3f007d")),
                                                  col=colorRamp2(quantile(x[[i]], c(0.00,0.25,0.50,0.75,0.9,0.99)),
                                                                 col=c("#c994c7","#fde0dd","#fcc5c0", "#fa9fb5", "#ae017e","#7a0177")), #for pol2
                                                  column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                                  show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE), ...
                                        ) ## legend will be shown only if the name is given for a heatmap.
                    }
          }

          return(ehml)
}

set.seed(123)
partition = paste0("", kmeans(xx$norm_matrix[[1]], centers = num_cluster)$cluster)
print(table(partition))
colr <- num_cluster+1
ehm_list <- get_enrichment_heatmap_list(x = xx_1$norm_matrix,names = xx_1$names,titles = xx_1$names,
                                        row_order = order(width(genes), decreasing = FALSE),
                                        axis_name_rot = 90,
                                        axis_name = c("-1kb","TSS", "+3kb"),
                                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col =2:colr), axis_param =list(facing="outside",side="left",
                                                                                                                  gp = gpar(fonsize=12, fontfamily="Arial")),
                                                                                                 ylim = c(0,20), height = unit(2.5, "cm")))
                                                                           )

ht_list = Heatmap(partition, col = structure(2:colr),
                  name="cluster",show_row_names = FALSE, width = unit(3, "mm"),
                  row_order = order(width(genes), decreasing = FALSE))+ehm_list


print("plotting....")
pdf(file=paste(output, length(genes), "hm.pdf", sep="_"), width=6, height=15)
draw(ht_list, split=partition,  heatmap_legend_side = "bottom", gap = unit(2.5, "mm"))
dev.off()








ehm_list_2 <- get_enrichment_heatmap_list(x = xx_2$norm_matrix,names = xx_2$names,titles = xx_2$names,
                                        row_order = order(width(genes), decreasing = FALSE),
                                        axis_name_rot = 90,
                                        axis_name = c("-1kb","TSS", "+3kb"),
                                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col =2:colr), axis_param =list(facing="outside",side="left",
                                                                                                                                          gp = gpar(fonsize=12, fontfamily="Arial")),
                                                                                                 ylim = c(0,40), height = unit(2.5, "cm")))
                                        )

ht_list_2 = Heatmap(partition, col = structure(2:colr),
                  name="cluster",show_row_names = FALSE, width = unit(3, "mm"),
                  row_order = order(width(genes), decreasing = FALSE))+ehm_list_2

print("plotting....")
pdf(file=paste("Pol2", length(genes), "hm.pdf", sep="_"), width=9, height=15)
draw(ht_list_2, split=partition,  heatmap_legend_side = "bottom", gap = unit(2.5, "mm"))
dev.off()

