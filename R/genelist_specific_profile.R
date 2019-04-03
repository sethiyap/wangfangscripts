#--- with TBP and Pol2

gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"
list_1 <- read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE) #-- plot sorted according to this list, provide expression value as well
list_2 <- read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE)
setwd(".") #make the folder with bw file as working directory
bw_test <- "H3K4me3_CGGACGTGG_CL_veA_wt_spore_ChIPmix22_normalized.bw" # specify the name of bw files
bw_control <- "H3_AGAACACC_an_spore_CL1019Mix_normalized.bw"
#--- for fig2D

#' genelist_specific_profiles
#'
#' plots two profiles for each provided genelist with normalising by the control
#' @param gff_file Provide reference gff file for your species of interest
#' @param bw_test Provide the bw file name to be tested
#' @param bw_control Provide bw file name to be normalised with i.e. H3
#' @param list_1 list of the genes of interest with expression value in second column
#' @param list_2 only one column of random or control genes
#' @param output heatmap in the pdf format
#'
#' @return
#' @export
#'
#' @examples
genelist_specific_profiles <- function(gff_file, bw_test,bw_control,list_1,list_2, output){

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

          #--- provide the genelist data with expression value for the list one

          list_1 <- list_1 %>% arrange(desc(X2))
          genes_1 <-  subset(genes, genes$gene_id %in% list_1$X1)
          genes_1 <- genes_1[match(list_1$X1,genes_1$gene_id),]

          #--- get second list co-ordinates
          genes_2 <-  subset(genes, genes$gene_id %in% list_2$X1)

          ## prepare signal data
          gene_lists <- list(genes_1, genes_2)

          names(gene_lists) <- paste(gsub(pattern = "_[[:upper:]]{6,}_.*_normalized.bw",replacement = "", bw_test),"_list",seq(1:length(gene_lists)), sep="")

          gene_lists <- tibble(name=names(gene_lists), data=gene_lists)
          print(gene_lists)

          bw_file_test <- import.bw(bw_test)

          bw_file_control=import.bw(bw_control)
          ## generate normalised matrix in tidy way


         dd <-  gene_lists %>%
                    dplyr::mutate(mat=purrr::map(data, function(i){
                    nn <- EnrichedHeatmap::normalizeToMatrix(bw_file_test, i, value_column = "score",background = 0,
                                                             smooth = TRUE,extend = c(1000))
                    nn[nn<0]=0
                    return(nn)
                    })) %>%

                    dplyr::mutate(mat_h3=purrr::map(data, function(i){
                             nn <- EnrichedHeatmap::normalizeToMatrix(bw_file_control, i, value_column = "score",background = 0,
                                                                      smooth = TRUE,extend = c(1000))
                             nn[nn<0]=0
                             return(nn)
                   }))

          #--- normalise by h3 ----
         dd2 <- dd  %>% mutate(norm_mat = map2(mat,mat_h3, function(x,y){ mm = (x+0.01) / (y+0.01)
                             return(mm)
                   } ))

          print(dd2)

          get_enrichment_heatmap_list <- function(x, names, titles, ...) {

                    ll <- length(x)

                    ## first heatmap
                    ehml <- EnrichedHeatmap(mat = x[[1]], name = names[[1]], column_title = titles[[1]], show_heatmap_legend = T,
                                            col=colorRamp2(quantile(x[[1]], c(0.1,0.5,0.6,0.9,0.99)),
                                                           col=c("#feebe2","#fcc5c0","#fa9fb5","#c51b8a","#7a0177")),
                                            use_raster = TRUE, ...)

                    ## several other heatmaps if length of x > 1.
                    if (ll > 1) {
                              for (i in 2:ll) {
                                        print(i)
                                        ehml <- ehml +
                                                  EnrichedHeatmap(
                                                            mat = x[[i]],
                                                            col=colorRamp2(quantile(x[[1]], c(0.1,0.5,0.6,0.9,0.99)),
                                                                           col=c("#feebe2","#fcc5c0","#fa9fb5","#c51b8a","#7a0177")),
                                                            name = ifelse(length(names) >= i, names[i], "NA"),
                                                            use_raster = TRUE,
                                                            column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                                            show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE), ...
                                                  ) ## legend will be shown only if the name is given for a heatmap.
                              }
                    }

                    return(ehml)
          }

          ehm_list <- get_enrichment_heatmap_list(x = dd2$norm_mat,names = dd2$name,
                                                  titles = dd2$name,
                                                  cluster_rows = FALSE,
                                                  row_order=NULL,
                                                  show_row_names = TRUE,
                                                  axis_name_rot = 90,
                                                  heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", legend_width = unit(3, "cm"),
                                                                              title_position = "topcenter",labels_gp = gpar(fonsize=12, fontfamily="Arial")),
                                                  axis_name = c("-1kb","TSS","TES", "+1kb"),
                                                  axis_name_gp = gpar(fonsize=12, fontfamily="Arial"),
                                                  top_annotation = HeatmapAnnotation(lines = anno_enriched(axis_param =list( facing="inside",side="left",gp=gpar(fonsize=12, fontfamily="Arial")),
                                                                                                           ylim = c(0.8,4),height = unit(2, "cm")
                                                                                                          #ylim = c(0.7,1.8),height = unit(2, "cm")
                                                                                                         # ylim = c(0.1,5.8),height = unit(2, "cm")
                                                                                                           )
                                                                                     )
          )

          # row_order_list = row_order(ehm_list)
          # list_1$H3K4me3 <- list_1[row_order_list,]$X1


          print("plotting....")
          pdf(file=paste(output, length(genes_1), "hm.pdf", sep="_"), width=6, height=60)
          draw(ehm_list, heatmap_legend_side = "top", gap = unit(2, "mm"))
          dev.off()

}









