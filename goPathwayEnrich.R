calculatePvalue <- function(data ,data_ ,Gene.Coverage ,adjust.method ) {
  
  data <- data[which(data$tad_name != "NA"),]
  tads <- unique(data$tad_name)
  
  data.with.p <- data.table(Term = character(),
                            TAD = character(),
                            P.value = numeric(),
                            stringsAsFactors = F)
  
  iterations <- c(1:length(tads))
  
  for (i in iterations){
    
    tad <- tads[i]
    tad.terms <- data[which(data$tad_name == tad),]
    terms.number <- dplyr::count(tad.terms,Term)
    iter <- c(1:nrow(terms.number))
    
    for (j in iter){
      
      hitInSample <- terms.number$n[j]
      hitInPop <- tad.terms[which(tad.terms$Term == as.character(terms.number[j,1])),denominator]
      failInPop <- Gene.Coverage -  hitInPop
      tad.genes <- data_[which(data_$tad_name == tad),]
      tad.genes <- tad.genes[which(tad.genes$Gene_id != ""),]
      sampleSize <- length(tad.genes$Gene_id)
      
      #test for over-representation, enrichment
      p.value <- phyper(hitInSample-1,hitInPop, failInPop,sampleSize, lower.tail=FALSE)
      
      temp <- data.table(Term = as.character(terms.number[j,1]),
                         TAD = as.character(tad),
                         P.value = as.numeric(p.value),
                         stringsAsFactors = F)
      
      data.with.p <- rbind(data.with.p ,temp)
      
    }
  }
  
  data.with.p <- data.with.p[which(data.with.p$P.value != "NA"),]
  data.with.p <- unique(data.with.p)
  
  #p adjust
  data.with.p$P.adjust <- p.adjust(data.with.p$P.value, method = adjust.method)
  
  return(data.with.p) 
}

enrichAll <- function(biodata, dbs){
  
  #preparing data for enrichR
  data_selected <- biodata %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  data_selected <- data_selected[which((data_selected$Gene_id != "NA") & (data_selected$Gene_id != "")), ]
  data_selected <- data_selected %>%
    dplyr::select(Gene_id,tad_name) %>%
    unique()
  dataForEnrich <-c(data_selected$Gene_id)
  
  #Enrichment with GO Molecular Function terms, GO Biological Process terms and KEGG Pathways
  #using enrichr interface to connect to EnrichR
  enriched <- enrichr(dataForEnrich, dbs)
  
  enriched_MF <- as.data.table(enriched[[dbs[1]]])
  enriched_MF <- subset(enriched_MF, P.value < 0.05)
  enriched_BP <- as.data.table(enriched[[dbs[2]]])
  enriched_BP <- subset(enriched_BP, P.value < 0.05)
  enriched_KEGG <- as.data.table(enriched[[dbs[3]]])
  enriched_KEGG <- subset(enriched_KEGG, P.value < 0.05)
  
  newList <- list(GO.MF = enriched_MF,GO.BP = enriched_BP,KEGG = enriched_KEGG,data.with.genes = data_selected)
  return(newList)
}


enrichPerTAD <- function(biodata,dbs){
  
  #preparing data for enrichR
  full.tads <- biodata %>%
    dplyr::select(tad_name,Gene_id) %>%
    as.data.table()
  
  data_ <- full.tads %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  data_ <- data_[which((data_$Gene_id != "NA") & (data_$Gene_id != "")), ]
  data_ <- unique(data_)
  
  unique.tads <- unique(data_$tad_name)
  
  #data.tables to store the enrichment results
  enriched_MF <- data.table(Term = character(),
                            Overlap = character(),
                            Genes = character(),
                            stringsAsFactors = F)
  
  enriched_BP <- data.table(Term = character(),
                            Overlap = character(),
                            Genes = character(),
                            stringsAsFactors = F)
  
  enriched_KEGG <- data.table(Term = character(),
                              Overlap = character(),
                              Genes = character(),
                              stringsAsFactors = F)
  
  #iterations <- c(1:10)
  iterations <- c(1:length(unique.tads))
  
  for (i in iterations){
    
    data.per.tad <- data_[which(data_$tad_name == unique.tads[i]),]
    dataForEnrich <-c(data.per.tad$Gene_id)
    dataForEnrich <- unique(dataForEnrich)
    
    #Enrichment with GO Molecular Function terms, GO Biological Process terms abd KEGG Pathways
    #using enrichr interface to connect to EnrichR
    enriched <- enrichr(dataForEnrich, dbs)
    
    loops <- c(1:length(dbs))
    for (l in loops){
      
      enriched_terms <- as.data.table(enriched[[dbs[l]]])
      
      if (nrow(enriched_terms)>0){
        enriched_terms <- subset(enriched_terms, P.value < 0.05)
        enriched_terms <- enriched_terms %>%
          dplyr::select(Term,Overlap,Genes)
        
        if(dbs[l] == "GO_Molecular_Function_2018"){
          enriched_MF <- rbind(enriched_MF,enriched_terms)
        }else if (dbs[l] == "GO_Biological_Process_2018"){
          enriched_BP <- rbind(enriched_BP,enriched_terms)
        }else{
          enriched_KEGG <- rbind(enriched_KEGG,enriched_terms)
        }
      }
    }
  }
  
  newList <- list(GO.MF = enriched_MF,GO.BP = enriched_BP,KEGG = enriched_KEGG,data.with.genes = data_)
  return(newList)
  
}

produceOutputs <- function(data.with.p,type){
  if ( str_detect(type,"GO")){
    
    data.with.p$Term <- str_remove(data.with.p$Term, "\\)")
    data.with.p <- data.with.p %>% 
      separate(Term, c("GO.term", "GO.number"), "\\(GO:" ) %>%
      as.data.table()
    
    data.visual <- data.with.p%>%
      dplyr::select(TAD, GO.term, GO.number, P.value, P.adjust)
    
    data.with.p<- data.visual %>%
      dplyr::group_by(TAD) %>%
      dplyr::summarise(TAD,GO.term = paste(GO.term, collapse = "|"),
                       GO.number = paste(GO.number, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    data.per.term <- data.visual %>%
      dplyr::group_by(GO.term, GO.number) %>%
      dplyr::summarise(GO.term, GO.number, 
                       TAD = paste(TAD, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    column.term <- paste0(type,".Term")
    column.id <- paste0(type,".number")
    column.p <- paste0(type,".P.value")
    column.adj <- paste0(type, ".P.adjust")
    colnames(data.with.p) <- c("TAD",column.term,column.id,column.p,column.adj)
    colnames(data.visual) <- c("TAD","Term","ID","P.value","P.adjust")
    colnames(data.per.term) <- c("GO.Term","GO.ID","TAD","P.value","P.adjust")
    
  }else{
    
    data.visual <- data.with.p%>%
      dplyr::select(TAD, Term, P.value, P.adjust)
    
    data.with.p<- data.visual %>%
      dplyr::group_by(TAD) %>%
      dplyr::summarise(TAD,Term = paste(Term, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    data.per.term <- data.visual %>%
      dplyr::group_by(Term) %>%
      dplyr::summarise(Term , 
                       TAD = paste(TAD, collapse = "|"),
                       p.value = paste(P.value, collapse = "|"),
                       P.adjust = paste(P.adjust, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    
    colnames(data.with.p) <- c("TAD","Term","P.value","P.adjust")
    colnames(data.visual) <- c("TAD","Term","P.value","P.adjust")
    colnames(data.per.term) <- c("TAD","Term","P.value","P.adjust")
    
  }
  
  
  newList <- list(data.visual = data.visual, data.perTerm = data.per.term, data.withP = data.with.p)
  return(newList)
}


analysisAll <- function(enriched_terms, type,data_selected, genes.coverage, p.adjust.method){
    
    if (nrow(enriched_terms)>0){
      
      enriched_terms <- enriched_terms %>%
        dplyr::select(Term,Overlap,Genes)
      
      #calculate number of genes per term in database 
      enriched_terms <- enriched_terms %>%
        separate(Overlap, c("numerator", "denominator"), sep = "\\/")
      enriched_terms$numerator <- as.numeric(as.character(enriched_terms$numerator))
      enriched_terms$denominator <- as.numeric(as.character(enriched_terms$denominator))
      
      data.extended <- enriched_terms %>%
        separate_rows(Genes, sep = ";", convert = TRUE) %>%
        as.data.table()
      
      data.extended <- left_join(data.extended,data_selected, by =c("Genes" = "Gene_id"))
      
      data.with.p <- data.table(Term = character(),
                                TAD = character(),
                                P.value = numeric(),
                                P.adjust = numeric(),
                                stringsAsFactors = F)
      
      data.with.p <- calculatePvalue(data = data.extended,data_ = data_selected, Gene.Coverage = genes.coverage, adjust.method = p.adjust.method)
      
      resultsList <- produceOutputs(data.with.p,type)
      
      return(resultsList)
    }
  
}


analysisPerTAD <- function(enriched_terms,type, data_, genes.coverage, p.adjust.method){

    enriched_terms <- enriched_terms %>%
      dplyr::select(Genes,Term,Overlap)
    
    #calculate number of genes per term in database
    enriched_terms <- enriched_terms %>%
      separate(Overlap, c("numerator", "denominator"), sep = "\\/")
    enriched_terms$numerator <- as.numeric(as.character(enriched_terms$numerator))
    enriched_terms$denominator <- as.numeric(as.character(enriched_terms$denominator))
    
    enriched_terms <- enriched_terms %>%
      dplyr::select(Term, denominator, Genes)%>%
      group_by(Term) %>%
      summarise(Term, 
                denominator = max(denominator),
                Genes = paste(Genes, collapse =";"),) %>%
      as.data.table() %>%
      unique()
    
    data.extended <- enriched_terms %>%
      separate_rows(Genes, sep = ";", convert = TRUE) %>%
      as.data.table()
    
    data.extended <- left_join(data.extended,data_, by =c("Genes" = "Gene_id")) %>%
      unique()
    data.extended <- data.extended %>%
      dplyr::select(Genes, Term, denominator, tad_name)
    
    data.with.p <- data.table(Term = character(),
                              TAD = character(),
                              P.value = numeric(),
                              P.adjust = numeric(),
                              stringsAsFactors = F)
    
    data.with.p <- calculatePvalue(data = data.extended,data_ = data_, Gene.Coverage = genes.coverage, adjust.method = p.adjust.method)
    
    resultsList <- produceOutputs(data.with.p,type)
    
    return(resultsList)
  
}

getKEGGIds <- function(enriched_KEGG){
  
  enriched_KEGG <- enriched_KEGG %>%
    dplyr::select(Term, Genes)
  enriched_KEGG <- enriched_KEGG %>% 
    group_by(Term) %>%
    summarise(Term,
              Genes = paste(Genes, collapse = ";"),) %>%
    as.data.table()
  enriched_KEGG <- unique(enriched_KEGG)
  
  #data for Pathview hsa ids
  data(paths.hsa)
  
  #get their kegg ids and save for later
  pathview.data <- data.table(hsa.ids = names(paths.hsa),
                              Term = paths.hsa,
                              stringsAsFactors = FALSE)
  pathview.data <- merge(pathview.data,enriched_KEGG,by = "Term") %>%
    dplyr::select(hsa.ids,Term,Genes)
  
  pathview.data <- unique(pathview.data)
  
  return(pathview.data)
}


createFolders <- function(output_folder){

  go_output_folder = paste0(output_folder,"/GO EA Outputs")
  kegg_output_folder = paste0(output_folder,"/Pathways EA Outputs")
  go_all_outputs = paste0(go_output_folder,"/GO enrich all")
  go_all_images = paste0(go_all_outputs,"/Images")
  go_per_outputs = paste0(go_output_folder,"/GO enrich per TAD")
  go_per_images = paste0(go_per_outputs,"/Images")
  kegg_all_outputs = paste0(kegg_output_folder,"/Path enrich all")
  kegg_all_images = paste0(kegg_all_outputs,"/Images")
  kegg_per_outputs = paste0(kegg_output_folder,"/Path enrich per TAD")
  kegg_per_images = paste0(kegg_per_outputs,"/Images")
  motif_output_folder = paste0(output_folder,"/Motif EA Outputs")
  image_output_folder = paste0(motif_output_folder,"/Images")
  
  dir.create(output_folder, showWarnings = FALSE)
  dir.create(go_output_folder, showWarnings = FALSE)
  dir.create(go_all_outputs, showWarnings = FALSE)
  dir.create(go_all_images, showWarnings = FALSE)
  dir.create(go_per_outputs, showWarnings = FALSE)
  dir.create(go_per_images, showWarnings = FALSE)
  dir.create(kegg_output_folder, showWarnings = FALSE)
  dir.create(kegg_all_outputs, showWarnings = FALSE)
  dir.create(kegg_all_images, showWarnings = FALSE)
  dir.create(kegg_per_outputs, showWarnings = FALSE)
  dir.create(kegg_per_images, showWarnings = FALSE)
  dir.create(motif_output_folder, showWarnings = FALSE)
  dir.create(image_output_folder, showWarnings = FALSE)
  
  folder.names <-c(go_all_outputs, go_all_images,go_per_outputs,go_per_images,
                   kegg_all_outputs,kegg_all_images,kegg_per_outputs,kegg_per_images,motif_output_folder,image_output_folder)  
  
  return(folder.names)
}

