

library(tidyverse)

dat <- tibble::tribble(
             ~Gene, ~RNA Pol-II Signal, ~mRNA_Expression,
            "atfA",            56.3248,      1051.890539,
              "wA",            36.6445,      826.2410956,
          "nce102",           216.3155,      819.7359104,
            "dewA",            21.5934,      453.8116876,
            "pkaA",            26.6596,      362.0308304,
            "pkaR",            25.7501,      340.6535569,
            "hogA",            50.1029,      336.5667278,
            "orlA",            21.6544,      253.5170659,
            "rhbA",            44.2994,      243.9687732,
            "wetA",            42.2475,      197.5830859,
            "vosA",            19.2343,      173.2160426,
          "AN5021",             28.038,      123.1145347,
          "AN6898",            19.1344,      97.24906504,
            "swoM",            14.5318,      95.30046085
          )


back2back_barplot <- function(dat, output){

          library(extrafont)
          loadfonts(device = "pdf")

          dat_m <- dat %>%
                    gather(Category, value, -Gene) %>%
                    mutate(value=log2(value), value=if_else(Category=="RNA Pol-II Signal", -value, value)) %>%
                    arrange(desc(value)) %>%
                    mutate(Gene=as_factor(Gene))

          gg <- ggplot(dat_m,aes(Gene, value, fill=Category))+geom_col()+
                    geom_hline(yintercept = 0)+
                    scale_fill_manual(values=c("blue1","#ef6548"))+
                    #scale_y_discrete(breaks=c(-10,-5,0,5,10),labels=c(10,5,0,5,10))+
                    coord_flip() +theme_bw()+ylim(c(-11,11))+
                    theme(axis.text.x= element_text(color="black",size=12,family="Arial",angle=0,vjust=0.8),
                          axis.text.y = element_text(color="black",size=12,family="Arial"),
                          axis.title.y=element_text(color="black",size=12,family="Arial", face="bold"),
                          axis.title.x=element_blank(),
                          legend.title=element_text(color="black",size=12,family="Arial"),
                          legend.key.size = unit(1,"line"),
                          legend.box.spacing=unit(0.5,"mm"),
                          legend.text=element_text(color="black",size=12,family="Arial"))+labs(y="log2(expression value)", x="")+
                    guides(fill=guide_legend(""))

          print(gg)
          ggsave(filename = paste(output,"_back2back_barplot.pdf",sep=""), device = "pdf", path="./", width = 18, height =15,unit="cm", dpi=300)




}
