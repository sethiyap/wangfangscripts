
dat <- clipr::read_clip()



#' plot_cnet
#'
#' @param dat provide vector of genes
#' @param output_name provide name of output file
#' @param org organism
#' @param GO_term_no no. of GO terms to be displayed
#' @import clusterprofiler
#' @import AnnotationHub
#' @import tidyverse
#' @return
#' @export
#'
#' @examples
plot_cnet <- function(dat, output_name, org,GO_term_no){

          library(tidyverse)
          library("clusterProfiler")
          library(AnnotationHub)
          #---
          dat <- unique(dat)
          #---
          data_provider <-"fungidb"  #"ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/"
          filter_tags = c(data_provider)  


          hub <- AnnotationHub()
          hub_subset <- query(hub , filter_tags)

          ## prepare summary table for annotation objects

          ah_data_summary <- tibble(ah_id = hub_subset$ah_id ,
                                    title = hub_subset$title ,
                                    dataprovider = hub_subset$dataprovider ,
                                    species = hub_subset$species ,
                                    taxonomyid = hub_subset$taxonomyid ,
                                    genome = hub_subset$genome ,
                                    description = hub_subset$description ,
                                    coordinate_1_based = hub_subset$coordinate_1_based ,
                                    maintainer = hub_subset$maintainer ,
                                    rdatadateadded = hub_subset$rdatadateadded ,
                                    preparerclass = hub_subset$preparerclass ,
                                    tags = hub_subset$tags ,
                                    rdataclass = hub_subset$rdataclass ,
                                    rdatapath = hub_subset$rdatapath,
                                    sourceurl = hub_subset$sourceurl ,
                                    sourcetype = hub_subset$sourcetype) #%>% filter(ah_id != "AH65059")

          ## re arrange summary data
          ah_data_summary <- ah_data_summary %>%
                    dplyr::select(ah_id, genome, species, taxonomyid,rdataclass) %>% spread(key = rdataclass, value = ah_id) %>% drop_na()

          species <- ah_data_summary %>% filter(species==org) %>% dplyr::select(OrgDb) %>% as.character()



          ### dowload orgdb
          my_species_orgdb <- hub_subset[[species]]

          print(my_species_orgdb)
          keys(my_species_orgdb)
          orgdb_cols <- c("GID", "GENEDESCRIPTION", "GO_ID" , "GO_TERM_NAME", "EVIDENCE_CODE" , "ONTOLOGY")

          my_species_orgdb_cols <- biomaRt::select(my_species_orgdb ,
                                                   keys = keys(my_species_orgdb) ,
                                                   columns = orgdb_cols, keytypes = "GID") %>% as_tibble()




          ## TERM to GENE table
          term_to_gene <- my_species_orgdb_cols  %>%
                    filter(ONTOLOGY == "Biological Process") %>%
                    dplyr::select(c("GO_ID", "GID")) %>% drop_na() %>% unique()

          ## TERM to NAME table
          term_to_name <- my_species_orgdb_cols %>% filter(ONTOLOGY == "Biological Process") %>%
                    dplyr::select(c("GO_ID", "GO_TERM_NAME")) %>% drop_na() %>% unique()

          enr <- enricher(gene = dat,
                          pvalueCutoff = 1,
                          pAdjustMethod = "none",
                          minGSSize = 1,
                          maxGSSize = 50000,
                          qvalueCutoff = 1,
                          TERM2GENE = term_to_gene,
                          TERM2NAME = term_to_name)

          enr@result %>% as_tibble() %>% clipr::write_clip()

          dotplot(enr, showCategory  = 10)

          print(emapplot(enr))

          cnetplot_out <- cnetplot(enr ,
                                   node_label = FALSE,
                                   colorEdge = TRUE ,
                                   showCategory = GO_term_no,
                                   layout = "auto")

          cnetplot_out <- cnetplot_out %+%
                    ggrepel::geom_text_repel(data = subset(cnetplot_out$data,color != "#B3B3B3") ,
                                             aes(x = x, y = y , label = name) ,
                                             size = 4) +
                    guides(edge_colour = "none")


          print(cnetplot_out)


          ggsave(filename = paste(output_name,"_clusterGO.pdf",sep=""), path="./", width = 18, height =10, dpi=300)



}
