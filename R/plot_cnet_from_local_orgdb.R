
#--- orgDB folder ----
dir="/Users/Pooja/Documents/Data-Analysis/PhD/Scripts-commands/orgDb/org.Scerevisiae.S288c.v42.eg.db_2019.03.tar.gz"
#--- Install packages for OrgDB ----
install.packages(pkgs = dir ,repos = NULL, type = "source")
#---


dat <- clipr::read_clip()

plot_cnet_from_local_orgdb <- function(dat, no_GO_Terms=20,output_name){

           library(tidyverse)
           library("clusterProfiler")
           library(AnnotationHub)

           dat <- unique(dat)

           #--- choose ordering by
           readinteger <- function()
           {
                     n <- readline(prompt="Enter index of species: ")
                     n <- as.integer(n)
                     if (is.na(n)){
                               n <- readinteger()
                     }
                     return(n)
           }


           packages <- installed.packages(.Library) %>% as_tibble() %>% dplyr::select(Package) %>% filter(str_detect(Package,"org"))

           print(packages)

           pkg_name <- (packages[readinteger(),])

           library(as.character(pkg_name),character.only = TRUE)

           #--- load org db
           my_species_orgdb <- eval(parse(text = pkg_name$Package))

           print(my_species_orgdb)

           #----- owrk on org db object ----
           orgdb_cols <- c("GID", "ANNOT_GENE_PRODUCT", "GO_GO_ID" , "GO_GO_TERM_NAME", "GO_EVIDENCE_CODE" , "GO_ONTOLOGY")

           my_species_orgdb_cols <- biomaRt::select(my_species_orgdb ,
                                                    keys = keys(my_species_orgdb) ,
                                                    columns = orgdb_cols, keytypes = "GID") %>% as_tibble()

           ## TERM to GENE table
           term_to_gene <- my_species_orgdb_cols  %>%
                     filter(GO_ONTOLOGY == "Biological Process") %>%
                     dplyr::select(c("GO_GO_ID", "GID")) %>% drop_na() %>% unique()

           ## TERM to NAME table
           term_to_name <- my_species_orgdb_cols %>% filter(GO_ONTOLOGY == "Biological Process") %>%
                     dplyr::select(c("GO_GO_ID", "GO_GO_TERM_NAME")) %>% drop_na() %>% unique()

           # term_to_gene <- term_to_gene %>% filter(GO_GO_ID!="GO:0008150") # remove  biological process GO
           # term_to_name <- term_to_name %>% filter(GO_GO_ID!="GO:0008150")

           #---- perform GO enrichment ----
           enr <- enricher(gene = dat,
                           pvalueCutoff = 1,
                           pAdjustMethod = "none",
                           minGSSize = 1,
                           maxGSSize = 50000,
                           qvalueCutoff = 1,
                           TERM2GENE = term_to_gene,
                           TERM2NAME = term_to_name)

           message("enrichment copied to clipboard")

           enr@result %>% as_tibble() %>% clipr::write_clip()

           cnetplot_out <- cnetplot(enr ,
                                    node_label = FALSE,
                                    colorEdge = TRUE ,
                                    showCategory = 20,
                                    layout = "auto")



           cnetplot_out <- cnetplot_out %+%
                     ggrepel::geom_text_repel(data = subset(cnetplot_out$data,color != "#B3B3B3") ,
                                              aes(x = x, y = y , label = name) ,
                                              size = 4) +
                     guides(edge_colour = "none")


           print(cnetplot_out)



           ggsave(filename = paste(output_name,"_clusterGO.pdf",sep=""), path="./", width = 12, height =10, dpi=300)


 }
