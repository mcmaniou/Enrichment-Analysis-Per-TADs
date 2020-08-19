########## Loading libraries ########## 
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(MotifDb)
library(ggseqlogo)
library(seqinr)
library(httr)
library(jsonlite)
library(xml2)
library("enrichR")
library(stats)
library(purrr)
library(igraph)
library(ggraph)
library(hrbrthemes)
library(extrafont)
library(RColorBrewer)
library(pathview)
library(gridExtra)

start_time = Sys.time()

source("motifEnrich.R")
source("goPathwayEnrich.R")
source("visualization.R")

#databases used from EnrichR
dbs<- c("GO_Molecular_Function_2018","GO_Biological_Process_2018","KEGG_2019_Human")

#gene coverage of the databases
genes.cover <- c(11459,14433,7802)

#choose a p adjust method, the methods supported are:
#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
p.adjust.method <- "fdr"

#RStudio supports different fonts for different operating systems
system = "win"

########### Inputs ########## 

dir_name = "Datasets"
output_folder = "Outputs"

folder.names <- create.folders(output_folder)

goAllOutputs <- folder.names[1]
goAllImages <- folder.names[2]
goPerOutputs <- folder.names[3]
goPerImages <- folder.names[4]
keggAllOutputs <- folder.names[5]
keggAllImages <- folder.names[6]
keggPerOutputs <- folder.names[7] 
keggPerImages <- folder.names[8]
motifOutputs <- folder.names[9]
motifImageOutputs <- folder.names[10]  

biodata = fread(paste(dir_name, "/integrated_table_with_sign_tads-ENSG.csv", sep = ""))

########### Enrichment + Data Analysis ##########

#enrichment all
listAll <- enrichAll(biodata,dbs)

#data analysis
data.type <- c("GO.MF","GO.BP","KEGG")

#GO MF Terms
data.MF <- data_analysis_all(listAll$GO.MF, data.type[1],listAll$data.with.genes, genes.cover[1], goAllOutputs,p.adjust.method)

#GO BP Terms
data.BP <- data_analysis_all(listAll$GO.BP, data.type[2],listAll$data.with.genes, genes.cover[2], goAllOutputs,p.adjust.method)

#KEGG Pathways
pathview.all <- getKEGGIds(listAll$KEGG)          #get pathview input data
data.KEGG.all <- data_analysis_all(listAll$KEGG, data.type[3],listAll$data.with.genes, genes.cover[3], keggAllOutputs,p.adjust.method)

#join GO Molecular Function and Biological Process outputs
data_all <- full_join(data.MF, data.BP, by = "TAD")


#enrichment per TAD
listPerTAD <- enrichPerTAD(biodata, dbs)

#GO MF Terms
data.MF <- analysis_perTAD(listPerTAD$GO.MF, data.type[1],listPerTAD$data.with.genes, genes.cover[1], goPerOutputs,p.adjust.method)

#GO BP Terms
data.BP <- analysis_perTAD(listPerTAD$GO.BP, data.type[2],listPerTAD$data.with.genes, genes.cover[2], goPerOutputs,p.adjust.method)

#KEGG Pathways
pathview.perTAD <- getKEGGIds(listPerTAD$KEGG)        #get pathview input data
data.KEGG.perTAD <- analysis_perTAD(listPerTAD$KEGG, data.type[3],listPerTAD$data.with.genes, genes.cover[3], keggPerOutputs,p.adjust.method)

#join GO Molecular Function and Biological Process outputs
data_perTAD <- full_join(data.MF, data.BP, by = "TAD")

#motif enrichment
listMotif <- motif_enrich(biodata, motifOutputs,p.adjust.method)

########### Output Files ########## 
end_enrich_all = Sys.time()

fwrite(data_all, paste(goAllOutputs, "/over-represented GO terms-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(data.KEGG.all, paste(keggAllOutputs, "/over-represented KEGG Pathways-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(data_perTAD, paste(goPerOutputs, "/over-represented GO terms-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(data.KEGG.perTAD, paste(keggPerOutputs, "/over-represented KEGG Pathways-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMotif$table_perTAD, paste0(motifOutputs, "/over-represented TFs in each tad.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMotif$table_perTFs, paste0(motifOutputs, "/TFs in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)

########### Visualization ##########

setGraphFonts(system)

#enrich all visualization
enrichrVisual(goAllOutputs , goAllImages)
enrichrVisual(keggAllOutputs , keggAllImages)
pathVisual(biodata, pathview.all ,keggAllImages)

#enrich per TAD visualization
enrichrVisual(goPerOutputs , goPerImages)
enrichrVisual(keggPerOutputs , keggPerImages)
pathVisual(biodata, pathview.perTAD ,keggPerImages)

#motif enrichment analysis visualization
motifVisual(motifImageOutputs, motifOutputs)

total_time <- Sys.time()


