mergeSequences <- function(data){
  
  #filter for the sequences related to TFs location
  data <- data[str_detect(data$Gene_locus,"promoter|intergenic|fiveUTR"),]
  data <- data %>%
    dplyr::select(tad_name,chromosome_name,start_position,end_position)
  
  #expand the CG sequences
  iter <- c(1: nrow(data))
  for (j in iter){
    if ((data$end_position[j] - data$start_position[j]) < 29){
      data$start_position[j] <- data$start_position[j] - 25
      data$end_position[j] <- data$end_position[j] + 25
    }
  }
  
  #group data according to TAD 
  TADs <- data %>%
    group_by(tad_name)
  groups <- group_split(TADs)  
  
  #merge overlapping sequences
  iter <- c(1: length(groups))
  for (j in iter){
    
    temp <- as.data.table(groups[[j]])
    temp <- temp[base::order(start_position, end_position)]
    if (nrow(temp)>1){
      iterations <- c(1:(nrow(temp)-1))
      
      for (k in iterations){
        if ( (temp$end_position[k] > temp$start_position[k+1])){
          
          temp$end_position[k] <- max(temp$end_position[k],temp$end_position[k+1])
          temp$start_position[k] <- min(temp$start_position[k],temp$start_position[k+1])
          temp$end_position[k+1] <- temp$end_position[k]
          temp$start_position[k+1] <-temp$start_position[k]
          
          loops <- c(1:k)
          for (l in loops){
            if (temp$start_position[l] == temp$start_position[k]){
              temp$end_position[l] <- temp$end_position[k]
              temp$start_position[l] <- temp$start_position[k]
            }
          }
        }
      }
    }
    temp <- unique(temp)
    groups[[j]] <- temp
  }
  
  data <- dplyr::bind_rows(groups)
  
  return(data)
}



getSequences <- function(input.data, outputs_folder){
  
  #query Ensembl Rest Api to get the sequences 
  #per TAD so as not to lose the TAD information
  new.TADs <- input.data %>%
    group_by(tad_name)
  new.groups <- group_split(new.TADs)
  
  iterations <- c(1:10)
  k <- 10+1
  #iterations <- c(1:length(new.groups))
  #k <- length(new.groups)+1
  seq.tad.number <- data.table(start = numeric(k),
                               end = numeric(k),
                               tad = character(k),
                               stringsAsFactors = F)
  seq.tad.number$start[1] <- 1
  
  for (i in iterations){
    
    data <- new.groups[[i]]
    iter <- c(1:nrow(data))  
    seq <- data.table(dna.seq = character(),
                      tad = character(),
                      stringsAsFactors = FALSE)
    
    for (j in iter){
      start <- data$start_position[j]
      end <- data$end_position[j]
      chr <- data$chromosome_name[j]
      server <- "http://rest.ensembl.org"
      ext <- paste("/sequence/region/human/",chr,":",start, "..", end, ":","1" ,"?", sep = "")
      r <- httr::GET(paste(server, ext, sep = ""), content_type("text/plain"))
      #stop_for_status(r)
      x <- data.table(dna.seq = httr::content(r),
                      tad = data$tad_name[j], 
                      stringsAsFactors = FALSE)
      seq <- rbind(seq,x)
    } 
    
    if (!is.null(seq)){
      
      write.fasta(sequences = as.list(seq$dna.seq) , names = seq$tad, file.out = paste0(outputs_folder,"/seq_perTADs.fasta"), open = "a")
      seq.tad.number$end[i] <- (nrow(seq) + seq.tad.number$start[i] -1)
      seq.tad.number$start[i+1] <- (seq.tad.number$end[i] +1)
      seq.tad.number$tad[i] <- data$tad_name[1]
      
    }
  }
  
  seq.tad.number <- seq.tad.number[which(seq.tad.number$end != 0),]
  return(seq.tad.number)
}

motif_outputs <- function(report.list, motif_output_folder){
  
  csv_perTAD <- data.table(TAD = character(),
                           TFs = character(),
                           Motifs = character(),
                           P.value = character())
  
  csv_perTF <- data.table(TFs = character(),
                          TAD = character(),
                          P.value = character(),
                          Motifs = character())
  
  iterations <- c(1: length(report.list))
  for (i in iterations){
    temp <- report.list[[i]]@d
    
    #exclude uncharacterized motifs
    temp <- temp[str_detect(temp$id,"UW.Motif.",negate = TRUE),]
    perTAD <- temp %>%
      dplyr::select(tad,target,id,p.value) %>%
      dplyr::summarise(TAD = tad,TFs = paste(target, collapse = "|"),
                       Motifs = paste(id, collapse = "|"),
                       P.value = paste(p.value, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    csv_perTAD <- rbind(csv_perTAD,perTAD)
    
    perTF <- temp %>%
      dplyr::select(target,tad,p.value,id) %>%
      group_by(target,tad) %>%
      dplyr::summarise(target, tad,
                       P.value = paste(p.value, collapse = "|"),
                       Motifs = paste(id, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    colnames(perTF) <- c("TFs", "TAD", "P.value", "Motifs")
    csv_perTF <- rbind(csv_perTF,perTF)
  }
  
  data.visual <- csv_perTF %>%
    dplyr::select(TFs,TAD,P.value)
  
  csv_perTF <- csv_perTF %>%
    group_by(TFs) %>%
    dplyr::summarise(TFs,TAD = paste(TAD,collapse ="|"),
                     P.value = paste(P.value, collapse = "|"),
                     Motifs = paste(Motifs, collapse = "|"),) %>%
    as.data.table()%>%
    unique()
  
  iterations <- c(1:nrow(csv_perTF))
  for(i in iterations){
    
    temp <- csv_perTF[i] %>%
      dplyr::select(TFs,Motifs)
    temp <- temp %>%
      separate_rows(Motifs, sep = "\\|") %>%
      as.data.table()
    temp <- unique(temp)
    
    temp <- temp %>%
      group_by(TFs) %>%
      dplyr::summarise(TFs,
                       Motifs = paste(Motifs, collapse = "|"),) %>%
      as.data.table()%>%
      unique()
    csv_perTF[i,4] <- temp[1,2]
  }
  
  table.TFs.Motifs <- csv_perTF %>%
    dplyr::select(TFs,Motifs)
  data.visual <- left_join(data.visual,table.TFs.Motifs)
  
  fwrite(data.visual, paste0(motif_output_folder, "/data for TFs visualization.csv"), 
                 row.names = FALSE, sep = ",", quote = FALSE)
  newList <- list(table_perTAD = csv_perTAD,table_perTFs = csv_perTF)
  return(newList)
}


motif_enrich <- function(motif.data, motif_output_folder, p.adjustment.method){
  
  #number of cores available for motif enrichment analysis
  #N <- 3
  
  #speed up execution
  #registerCoresPWMEnrich(N)
  #useBigMemoryPWMEnrich(TRUE)
  
  start_motifEA_time <- Sys.time()
  
  motif.data <- mergeSequences(motif.data)
  
  merge_time = Sys.time()
  
  seq.tad.number <- getSequences(motif.data,motif_output_folder)
  
  translate_time <- Sys.time()
  
  #perform motif enrichment analysis using PWMEnrich
  report.list <- list()
  
  # load the pre-compiled lognormal background
  data(PWMLogn.hg19.MotifDb.Hsap)
  l <- 1
  iterations <- c(1:10)
  #iterations <- c(1:nrow(seq.tad.number))
  for (i in iterations){
    
    sequence = readDNAStringSet(paste0(motif_output_folder,"/seq_perTADs.fasta"), format="fasta", skip = (seq.tad.number$start[i]-1), nrec = (seq.tad.number$end[i]-seq.tad.number$start[i]+1) )
    res = motifEnrichment(sequence, PWMLogn.hg19.MotifDb.Hsap)
    report = groupReport(res, by.top.motifs = TRUE)
    report <- report[report$p.value < 0.05]
    report.list[[l]]<-report
    report.list[[l]]@d$tad <- seq.tad.number$tad[i]
    
    l <- l+1
    
  }
  
  end_time = Sys.time()
  
  #registerCoresPWMEnrich(NULL)
  #useBigMemoryPWMEnrich(FALSE)
  
  file.create(paste0(motif_output_folder,"/report_motif.txt"), showWarnings = FALSE)
  dput(paste0(motif_output_folder,"/report_motif.txt"), file = "report_motif.txt")
  
  #report.list <- dget(paste0(motif_output_folder,"/report_motif.txt"))
  
  listMotifEA <- motif_outputs(report.list, motif_output_folder)
  
  return(listMotifEA)
  
}