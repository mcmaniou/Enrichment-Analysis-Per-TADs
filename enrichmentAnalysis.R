########## Loading libraries ########## 
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(tidyverse)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggseqlogo)
library(seqinr)
library(httr)
library(jsonlite)
library(xml2)
library(enrichR)
library(stats)
library(purrr)
library(igraph)
library(ggraph)
library(hrbrthemes)
library(extrafont)
library(pathview)
library(gridExtra)
library(saveImageHigh)
library(KEGGREST)
library(ggpubr)

start_time = Sys.time()

source("motifEnrich.R")
source("goPathwayEnrich.R")
source("visualization.R")

#databases used from EnrichR
dbs <- c("GO_Molecular_Function_2018","GO_Biological_Process_2018","KEGG_2019_Human")

#gene coverage of the databases
genes.cover <- c(11459,14433,7802)

#choose a p adjust method, the methods supported are:
#c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
p.adjust.method <- "fdr"

#cut-off enrichment adjusted p-value
cut.off <- 0.05

#min number of genes in over-represented terms
min.genes <- 3

#RStudio supports different fonts for different operating systems
system = "win"

########### Inputs ########## 

dir_name = "Datasets"
output_folder = "Outputs"

folder.names <- createFolders(output_folder)

goAllOutputs <- folder.names[1]
goAllImages <- folder.names[2]
goPerOutputs <- folder.names[3]
goPerImages <- folder.names[4]
keggAllOutputs <- folder.names[5]
keggAllImages <- folder.names[6]
keggPerOutputs <- folder.names[7] 
keggPerImages <- folder.names[8]
motifOutputsFolder <- folder.names[9]
motifImageOutputs <- folder.names[10]  

biodata = fread(paste(dir_name, "/integrated_table_with_sign_tads-sample_input.csv", sep = ""))

########### Enrichment + Data Analysis ##########

#enrichment all
listAll <- enrichAll(biodata,dbs, cut.off)

#data analysis
data.type <- c("GO.MF","GO.BP","KEGG")

#GO MF Terms
listMFAll <- analysisAll(listAll$GO.MF, data.type[1],listAll$data.with.genes, genes.cover[1],p.adjust.method, min.genes)

#GO BP Terms
listBPAll <- analysisAll(listAll$GO.BP, data.type[2],listAll$data.with.genes, genes.cover[2],p.adjust.method, min.genes)

#KEGG Pathways
pathviewAll <- getKEGGIds(listAll$KEGG)          #get pathview input data
listKEGGAll <- analysisAll(listAll$KEGG, data.type[3],listAll$data.with.genes, genes.cover[3],p.adjust.method , min.genes)

#join GO Molecular Function and Biological Process outputs
dataAll <- full_join(listMFAll$data.withP, listBPAll$data.withP, by = "TAD")


#enrichment per TAD
listPerTAD <- enrichPerTAD(biodata, dbs, cut.off)

#GO MF Terms
listMFPerTAD <- analysisPerTAD(listPerTAD$GO.MF, data.type[1],listPerTAD$data.with.genes, genes.cover[1], p.adjust.method, min.genes)

#GO BP Terms
listBPPerTAD <- analysisPerTAD(listPerTAD$GO.BP, data.type[2],listPerTAD$data.with.genes, genes.cover[2], p.adjust.method, min.genes)

#KEGG Pathways
pathviewPerTAD <- getKEGGIds(listPerTAD$KEGG)        #get pathview input data
listKEGGPerTAD <- analysisPerTAD(listPerTAD$KEGG, data.type[3],listPerTAD$data.with.genes, genes.cover[3], p.adjust.method, min.genes)

#join GO Molecular Function and Biological Process outputs
dataPerTAD <- full_join(listMFPerTAD$data.withP, listBPPerTAD$data.withP, by = "TAD")

#motif enrichment
listMotif <- motifEnrich(biodata, motifOutputsFolder)

########### Output Files ########## 
end_enrich_all = Sys.time()

#enrichment all
fwrite(dataAll, paste(goAllOutputs, "/over-represented GO terms-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMFAll$data.perTerm, paste0(goAllOutputs, "/GO MF Terms in different TADs.csv"), 
      row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listBPAll$data.perTerm, paste0(goAllOutputs, "/GO BP Terms in different TADs.csv"), 
      row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGAll$data.withP, paste(keggAllOutputs, "/over-represented KEGG Pathways-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGAll$data.perTerm, paste0(keggAllOutputs, "/KEGG Pathways in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)

#enrichment per TAD
fwrite(dataPerTAD, paste(goPerOutputs, "/over-represented GO terms-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMFPerTAD$data.perTerm, paste0(goPerOutputs, "/GO MF Terms in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listBPPerTAD$data.perTerm, paste0(goPerOutputs, "/GO BP Terms in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGPerTAD$data.withP, paste(keggPerOutputs, "/over-represented KEGG Pathways-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGPerTAD$data.perTerm, paste0(keggPerOutputs, "/KEGG Pathways in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)

#motif enrichment
fwrite(listMotif$table_perTAD, paste0(motifOutputsFolder, "/over-represented TFs in each tad.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMotif$table_perTFs, paste0(motifOutputsFolder, "/TFs in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)

########### Visualization ##########

setGraphFonts(system)

#enrich all visualization
enrichrVisual(goAllImages,"GO MF Terms",listMFAll$data.visual)
enrichrVisual(goAllImages, "GO BP Terms",listBPAll$data.visual)
enrichrVisual(keggAllImages, "KEGG Pathways",listKEGGAll$data.visual)
pathVisual(biodata, pathviewAll ,keggAllImages, keggAllOutputs)

#enrich per TAD visualization
enrichrVisual(goPerImages,"GO MF Terms",listMFPerTAD$data.visual)
enrichrVisual(goPerImages,"GO BP Terms",listBPPerTAD$data.visual)
enrichrVisual(keggPerImages, "KEGG Pathways",listKEGGPerTAD$data.visual)
pathVisual(biodata, pathviewPerTAD ,keggPerImages, keggPerOutputs)

#motif enrichment analysis visualization
report.list <- dget(paste0(motifOutputsFolder,"/report_motif.txt"))
motifVisual(motifImageOutputs, motifOutputsFolder, listMotif$data.visual, report.list)

total_time <- Sys.time()


