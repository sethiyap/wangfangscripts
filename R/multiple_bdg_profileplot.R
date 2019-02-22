#--- with TBP and Pol2

gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"
mylist <- read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE)
bw_files_dir <-c(".")
#---

#' multiple_bdg_profileplot
#'
#' @param gff_file  Provide reference gff file for your species of interest
#' @param bw_files_dir Put all the bw files in owrking directory with ending like "_CACAGTTGG_CL1019Mix_normalized.bw"
#' @param output heatmap in the pdf format
#' @param mylist provide customised list if want to plot specific genes from the entire genome, else specify mylist="NULL"
#' @param cutoff provide cut-off number (ie. row number) to highlight region like actively transcribing genes, or specific gene if nothing to specific use 0
#'
#' @return
#' @export
#'
#' @examples
#'
multiple_bdg_profileplot <- function(gff_file, bw_files_dir, output, mylist, cutoff){

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

                    row_order = order(enriched_score(mat), decreasing = TRUE)
          }

          else{
                    genes <-  subset(genes, genes$gene_id %in% mylist$X1)

                    genes <- genes[match(mylist$X1,genes$gene_id),]

                    row_order = order(mylist$X2, decreasing=TRUE)
          }

          tss <-  promoters(genes,upstream=0, downstream=1)

          ## prepare signal data
          bw_files <- list.files(bw_files_dir, pattern = "*_normalized.*", recursive = T, full.names = T)

          #--- replace file names to get the sample names
          #eg: ./H2A.Z_by_Htz1_NH4_GTTGTCCCA_CW466_525_normalized.bw
          # print only H2A.Z_by_Htz1_NH4

          names(bw_files) <- gsub(pattern = "_[[:upper:]]{6,}_.*_normalized.bw",replacement = "", bw_files)
          names(bw_files) <- gsub(pattern = "./",replacement = "", names(bw_files))

          #names(bw_files) <- c("hepA_HA_20hr","hepA_HA_48hr","hypha_h3k9me3_wt","hypha_h3k9me3_hepA_del", "spore_h3k9me3_wt","spore_h3k9me3_hepA_del" )
          bw_files <- broom::tidy(bw_files)
          print(bw_files)

          #--- choose ordering by
          readinteger <- function()
          {
                    n <- readline(prompt="Enter index of file by which heatmap to be clustered: ")
                    n <- as.integer(n)
                    if (is.na(n)){
                              n <- readinteger()
                    }
                    return(n)
          }


          ## generate normalised matrix in tidy way
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


          get_enrichment_heatmap_list <- function(x, names, titles, ...) {
                    #x <- xx$norm_matrix

                    ll <- length(x)

                    ## first heatmap
                    ehml <- EnrichedHeatmap(mat = x[[1]], name = names[[1]], column_title = titles[[1]], show_heatmap_legend = T,
                                            use_raster = TRUE, ...)

                    ## several other heatmaps if length of x > 1.
                    if (ll > 1) {
                              for (i in 2:ll) {
                                        print(i)
                                        ehml <- ehml +
                                                  EnrichedHeatmap(
                                                            mat = x[[i]],
                                                            name = ifelse(length(names) >= i, names[i], "NA"),
                                                            use_raster = TRUE,
                                                            right_annotation = ha,
                                                            column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                                            show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE), ...
                                                  ) ## legend will be shown only if the name is given for a heatmap.
                              }
                    }

                    return(ehml)
          }



          ha = rowAnnotation(foo = anno_mark(at = cutoff, labels = "actively transcribing"))


          ehm_list <- get_enrichment_heatmap_list(x = xx$norm_matrix,names = xx$names,titles = xx$names,
                                                  cluster_rows = FALSE,
                                                  row_order = row_order,
                                                  show_row_names = FALSE,
                                                  axis_name_rot = 90,
                                                  heatmap_legend_param = list(color_bar = "continuous",legend_direction="vertical"),
                                                  axis_name = c("-1kb","TSS", "+1kb"),
                                                  axis_name_gp = gpar(fonsize=12, fontfamily="Arial"),
                                                  col = colorRamp2(breaks = seq(0, 60, by = 10),
                                                                   colors = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#b10026")),
                                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(axis_param =list( facing="inside",gp=gpar(fonsize=12, fontfamily="Arial"))))
                                                 )

          print("plotting....")
          pdf(file=paste(output, length(genes), "hm.pdf", sep="_"), width=6, height=10)
          draw(ehm_list, heatmap_legend_side = "right", gap = unit(2.1, "mm"))
          dev.off()

}
