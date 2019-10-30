
gff_file <- "/Users/Pooja/Documents/Data-Analysis/Others/Reference_Annotation/An/A_nidulans_FGSC_A4_version_s10-m04-r07_features.gff"
mylist <- readr::read_delim(pipe("pbpaste"), delim="\t", col_names =FALSE)
bw_files_dir <-c(".")
#---

#' multiple bdg profileplot
#' @description obtain profile plot with user defined cluster and ylim for the line-plot
#'
#' @param gff_file a string, defining absolute path of the gff file
#' @param bw_files_dir a string, defining absolute path of the bw files
#' @param output a string, defining output file name
#' @param num_cluster numeric, number of clusters to be displayed in the plot
#' @param mylist use defined list of genes, if subset of genes are to be plotted, default: NULL
#' @param ylim numeric,define the y-axis limit for lineplot
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures promoters
#' @importFrom GenomicFeatures genes
#' @importFrom rtracklayer import
#' @importFrom broom tidy
#' @importFrom dplyr mutate
#' @import ComplexHeatmap
#' @import EnrichedHeatmap
#' @import grid
#' @importFrom BiocGenerics width
#' 
#' @return
#' @export
#'
#' @examples
#' multiple_bdg_profileplot(gff_file =gff_file,bw_files_dir = ".",num_cluster = 1,mylist =mylist,output = "sample", ylim = 150)
#' 
#' 
multiple_bdg_profileplot <- function(gff_file, bw_files_dir, output, num_cluster, mylist=NULL, ylim=5){
          
          ##--- packages
           library(dplyr)
          
          gff <- GenomicFeatures::makeTxDbFromGFF(gff_file, metadata = T)
          genes <- GenomicFeatures::genes(gff)
          
          if(is.null(mylist)== TRUE){
                    genes <- genes
          }
          
          else{
               genes <-  subset(genes, genes$gene_id %in% mylist$X1)     
          }
          
          tss <-  GenomicFeatures::promoters(genes,upstream=0, downstream=1)
          
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
                    dplyr::mutate(bw = purrr::map(x, function(ii) {
                              rtracklayer::import(ii)
                    })) %>%
                    dplyr::mutate(norm_matrix = purrr::map(bw, function(ii) {
                              nn <- EnrichedHeatmap::normalizeToMatrix(ii, tss, value_column = "score", mean_mode = "w0",
                                                                       keep = c(0.0, 0.99),background = 0, 
                                                                       smooth = TRUE,extend = c(1000, 3000), w = 10)
                              nn[nn<0]=0
                              return(nn)
                    }))
          
          
          get_enrichment_heatmap_list <- function(x, names, titles, ...) {
                    #x <- xx$norm_matrix
                    
                    ll <- length(x)
                    
                    ## first heatmap
                    ehml <- EnrichedHeatmap::EnrichedHeatmap(mat = x[[1]], name = names[[1]], column_title = titles[[1]], show_heatmap_legend = T,
                                            use_raster = TRUE,
                                            col=circlize::colorRamp2(quantile(x[[1]], c(0.00,0.25,0.50,0.75,0.9,0.99)), 
                                                           col=c("turquoise1","gold","goldenrod3", "skyblue","slateblue", "black")), ...)
                    
                    ## several other heatmaps if length of x > 1.
                    if (ll > 1) {
                              for (i in 2:ll) {
                                        print(i)
                                        ehml <- ehml +
                                                  EnrichedHeatmap::EnrichedHeatmap(
                                                            mat = x[[i]],
                                                            col=circlize::colorRamp2(quantile(x[[i]], c(0.00,0.25,0.50,0.75,0.9,0.99)), 
                                                                           col=c("turquoise1","gold","goldenrod3", "skyblue","slateblue", "black")),
                                                            name = ifelse(length(names) >= i, names[i], "NA"),
                                                            use_raster = TRUE,
                                                            column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                                            show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE), ... 
                                                  ) ## legend will be shown only if the name is given for a heatmap.
                              }
                    }
                    
                    return(ehml)
          }
          
          set.seed(123)
          partition = paste0("", kmeans(xx$norm_matrix[[readinteger()]], centers = num_cluster)$cluster)
          print(table(partition))
          colr <- num_cluster+1
          ehm_list <- get_enrichment_heatmap_list(x = xx$norm_matrix,
                                                  names = xx$names,
                                                  titles = xx$names, 
                                                  row_order = order(BiocGenerics::width(genes), decreasing = FALSE),
                                                  axis_name_rot = 90,
                                                  heatmap_legend_param = list(color_bar = "continuous",legend_direction="horizontal", legend_width = grid::unit(3, "cm"),
                                                                              title_position = "topcenter",labels_gp = grid::gpar(fonsize=12)),
                                                  axis_name_gp = grid::gpar(fonsize=12),
                                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(lines = EnrichedHeatmap::anno_enriched(axis_param =list( facing="inside",side="left",gp=grid::gpar(fonsize=12)),
                                                                                                           ylim = c(0.1,ylim),
                                                                                                           height = grid::unit(2, "cm")
                                                  )
                                                  ))
          
          
         
                                                
          
          ht_list = ComplexHeatmap::Heatmap(partition, col = structure(2:colr),
                            name="cluster",show_row_names = FALSE, width = grid::unit(3, "mm"),  
                            row_order = order(BiocGenerics::width(genes), decreasing = FALSE))+ehm_list

           


           print("plotting....")
           pdf(file=paste(output, length(genes), "hm.pdf", sep="_"), width=16, height=15)
           ComplexHeatmap::draw(ht_list, split=partition,  heatmap_legend_side = "bottom", gap = grid::unit(1.5, "mm"))
           dev.off()
          
           print("writing clusters ......")
           genes$partition = partition
           write.table(genes, file=paste(output, length(genes), "cluster.tab", sep="_"), sep="\t", row.names = FALSE, quote = FALSE)
          
}
