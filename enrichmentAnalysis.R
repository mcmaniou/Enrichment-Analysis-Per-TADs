########## Loading libraries ########## 
library(pryr)
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

source("motifEnrich.R")
source("goPathwayEnrich.R")
source("visualization.R")

########### Inputs ########## 

#scal_test <- data.table(part = character(6),
#                        time = character(6),
#                        memory = numeric(6) )

#scal_test$part[1] <- "start"
#scal_test$time[1] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
#scal_test$memory[1] <- mem_used()

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

#dir_name = "Datasets"
#output_folder = "Outputs"
#filepath = paste(dir_name, "/integrated_table_with_sign_tads-sample_input.csv", sep = "")

dir_name = "Datasets-all"
output_folder = "Outputs-1"
filepath = paste(dir_name, "/integrated_table_with_sign_tads-ENSG.csv", sep = "")

folder <- createFolders(output_folder)

biodata = fread(filepath)

#biodata <- biodata[40000:63426,]

########### Enrichment + Data Analysis ##########
#scal_test$part[2] <- "before enrich all"
#scal_test$time[2] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
#scal_test$memory[2] <- mem_used()

#enrichment all
listAll <- enrichAll(biodata,dbs, cut.off)

#data analysis
data.type <- c("GO.MF","GO.BP","KEGG")

#GO MF Terms
listMFAll <- dataAnalysis(listAll$GO.MF, data.type[1],listAll$data.with.genes, genes.cover[1],p.adjust.method, min.genes)

#GO BP Terms
listBPAll <- dataAnalysis(listAll$GO.BP, data.type[2],listAll$data.with.genes, genes.cover[2],p.adjust.method, min.genes)

#KEGG Pathways
listPathAll <- getKEGGIds(listAll$KEGG, biodata[,Gene_id,diff])          #get pathview input data
listKEGGAll <- dataAnalysis(listAll$KEGG, data.type[3],listAll$data.with.genes, genes.cover[3],p.adjust.method , min.genes)

#join GO Molecular Function and Biological Process outputs
dataAll <- full_join(listMFAll$data.perTAD, listBPAll$data.perTAD, by = "TAD")

#scal_test$part[3] <- "before enrich per TAD"
#scal_test$time[3] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
#scal_test$memory[3] <- mem_used()
#enrichment per TAD
listPerTAD <- enrichPerTAD(biodata, dbs, cut.off)

#GO MF Terms
listMFPerTAD <- dataAnalysis(listPerTAD$GO.MF, data.type[1],listPerTAD$data.with.genes, genes.cover[1], p.adjust.method, min.genes)

#GO BP Terms
listBPPerTAD <- dataAnalysis(listPerTAD$GO.BP, data.type[2],listPerTAD$data.with.genes, genes.cover[2], p.adjust.method, min.genes)

#KEGG Pathways
listPathPerTAD <- getKEGGIds(listPerTAD$KEGG, biodata[,Gene_id,diff])        #get pathview input data
listKEGGPerTAD <- dataAnalysis(listPerTAD$KEGG, data.type[3],listPerTAD$data.with.genes, genes.cover[3], p.adjust.method, min.genes)

#join GO Molecular Function and Biological Process outputs
dataPerTAD <- full_join(listMFPerTAD$data.perTAD, listBPPerTAD$data.perTAD, by = "TAD")

#scal_test$part[4] <- "before motif EA"
#scal_test$time[4] <- format(Sys.time(), "%d-%b-%Y %H.%M.%S")
#scal_test$memory[4] <- mem_used()

#motif enrichment
report.list <- motifEnrich(biodata, folder$motifOutputsFolder,p.adjust.method, cut.off)
#report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
listMotif <- motifOutputs(report.list)

########### Output Files ########## 
#scal_test$part[5] <- "before outputs"
#scal_test$time[5] <- format(Sys.time(), "%d-%b-%Y %H.%M")
#scal_test$memory[5] <- mem_used()

#enrichment all
fwrite(dataAll, paste(folder$goAllOutputs, "/over-represented GO terms-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMFAll$data.perTerm, paste0(folder$goAllOutputs, "/GO MF Terms in different TADs.csv"), 
      row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listBPAll$data.perTerm, paste0(folder$goAllOutputs, "/GO BP Terms in different TADs.csv"), 
      row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGAll$data.perTAD, paste(folder$keggAllOutputs, "/over-represented KEGG Pathways-enrichment all.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGAll$data.perTerm, paste0(folder$keggAllOutputs, "/KEGG Pathways in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listPathAll$output.csv, paste0(folder$keggAllOutputs, "/Pathview input.csv"),
       row.names = FALSE, sep = "\t", quote = FALSE)


#enrichment per TAD
fwrite(dataPerTAD, paste(folder$goPerOutputs, "/over-represented GO terms-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMFPerTAD$data.perTerm, paste0(folder$goPerOutputs, "/GO MF Terms in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listBPPerTAD$data.perTerm, paste0(folder$goPerOutputs, "/GO BP Terms in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGPerTAD$data.perTAD, paste(folder$keggPerOutputs, "/over-represented KEGG Pathways-enrichment per tad.csv", sep = ""), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listKEGGPerTAD$data.perTerm, paste0(folder$keggPerOutputs, "/KEGG Pathways in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listPathPerTAD$output.csv, paste0(folder$keggPerOutputs, "/Pathview input.csv"),
       row.names = FALSE, sep = "\t", quote = FALSE)

#motif enrichment
fwrite(listMotif$table_perTAD, paste0(folder$motifOutputsFolder, "/over-represented TFs in each tad.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
fwrite(listMotif$table_perTFs, paste0(folder$motifOutputsFolder, "/TFs in different TADs.csv"), 
       row.names = FALSE, sep = "\t", quote = FALSE)
file.create(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"), showWarnings = FALSE)
dput(report.list, file = paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))

########### Visualization ##########

setGraphFonts(system)

#enrich all visualization
enrichrVisual(folder$goAllImages,"GO MF Terms",listMFAll$data.visual)
enrichrVisual(folder$goAllImages, "GO BP Terms",listBPAll$data.visual)
enrichrVisual(folder$keggAllImages, "KEGG Pathways",listKEGGAll$data.visual)
pathVisual(listPathAll$pathview.input ,folder$keggAllImages)

#enrich per TAD visualization
enrichrVisual(folder$goPerImages,"GO MF Terms",listMFPerTAD$data.visual)
enrichrVisual(folder$goPerImages,"GO BP Terms",listBPPerTAD$data.visual)
enrichrVisual(folder$keggPerImages, "KEGG Pathways",listKEGGPerTAD$data.visual)
pathVisual(listPathPerTAD$pathview.input ,folder$keggPerImages)

#motif enrichment analysis visualization
#report.list <- dget(paste0(folder$motifOutputsFolder,"/report MotifEA.txt"))
motifVisual(folder$motifImageOutputs, folder$motifOutputsFolder, listMotif$data.visual, report.list)

#scal_test$part[6] <- "end"
#scal_test$time[6] <- format(Sys.time(), "%d-%b-%Y %H.%M")
#scal_test$memory[6] <- mem_used()

#fwrite(scal_test, paste0(getwd(), "/scalability test 4-new.csv"),
#       row.names = FALSE, sep = "\t", quote = FALSE)
