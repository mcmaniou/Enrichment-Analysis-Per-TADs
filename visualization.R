setGraphFonts <- function(deviceSystem){
  
  if (deviceSystem == "win"){
    #font_import() 
    loadfonts(device = "win")
    windowsFonts(Times=windowsFont("TT Times New Roman"))
    theme_set(theme_bw(base_size=12, base_family = 'Times New Roman')+ 
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()))
  }
}


enrichrVisual <- function(input_folder, output_folder){
  
  csvs <- c("/data for GO.MF visualization.csv","/data for GO.BP visualization.csv", "/data for KEGG visualization.csv")
  dir <- c("GO MF Terms", "GO BP Terms", "KEGG Pathways")
  
  if (str_detect(input_folder,"GO")){
    loops <- c(1:2) 
  }else{
    loops <- 3
  }
  
  for (l in loops){
    
    data.visual <-fread(paste(input_folder,csvs[l], sep = ""))
    
    if (str_detect(dir[l],"GO")){
      string <- paste(output_folder,"/", dir[l], sep = "")
      dir.create(string, showWarnings = FALSE) 
    }else{
      string <- output_folder
    }
    
    # plot 1
    data.plot1 <- data.table(P.value = c(data.visual$P.value, data.visual$P.adjust),
                             Type = c(rep_len("P value",nrow(data.visual)),rep_len("Adjusted P value",nrow(data.visual))),
                             stringsAsFactors = F)
    
    png(filename = paste(string, "/","Density plot-P values of ",dir[l], ".png", sep = ""), 
        width = 500, height = 650)
    
    p1 <- ggplot(data=data.plot1, aes(x=P.value, group=Type, fill=Type)) +
      geom_density(adjust=1.5, alpha=.4) +
      ggtitle(paste0("P values of " ,dir[l])) +
      xlab("P value")+
      ylab("Density")+
      theme_ipsum()+
      theme(
        plot.title = element_text(size=15),
        axis.title.x = element_text(size = 15, vjust = 0.5,hjust = 0.5),
        axis.title.y = element_text(size = 15, vjust = 1.5,hjust = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_text(size=15)
      )
    
    print(p1)
    dev.off()
    
    
    #plot 2
    terms <- dplyr::count(data.visual,Term)
    
    data.plot2 <- merge(terms, data.visual) %>%
      dplyr::select(Term, n, P.adjust) %>%
      group_by(Term,n) %>%
      summarise(Term, n, P.adjust = min(P.adjust)) %>%
      unique() 
    
    data.plot2 <- setorder(data.plot2,-n)
    
    if (nrow(data.plot2)>29){
      data.plot2 <- data.plot2[1:30,]
    }
    
    png(filename = paste(string, "/",dir[l]," in different TADs", ".png", sep = ""), 
        width = 900, height = 660)
    
    p2 <- data.plot2 %>%
      ggplot(aes(x = as.factor(reorder(Term, n)),y = n, fill = P.adjust)) +
      geom_histogram(color = "#e9ecef",stat = "identity") +
      labs(title = paste0("Top 30 ",dir[l]," in different TADs"), 
           subtitle = "Bar color corresponds to min Adjusted P value of each Term in the different TADs")+
      xlab(dir[l])+
      ylab("Number of TADs")+
      coord_flip() +
      theme_ipsum() +
      theme(
        plot.title = element_text(size=15),
        plot.subtitle = element_text(size=10),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12, vjust = 0.5,hjust = 0.5),
        axis.title.y = element_text(size = 12, vjust = 0.5,hjust = 0.5)
      )
    print(p2)
    dev.off()
    
    
    #plot 3
    #plot network graph only for GO terms
    if (str_detect(dir[l],"GO")){
     
      #choose only 10 largest groups
      go_count <- dplyr::count(data.visual,Term)
      go_count <- slice_max(go_count, order_by = n,n =10,with_ties = FALSE)
      data.plot5 <- merge(data.visual,go_count, by = "Term")
      
      nodes <- data.plot5%>%
        dplyr::select(TAD,Term)
      
      gos <- nodes %>% 
        group_by(Term)
      groups <- group_split(gos)
      
      iterations <- c(1:length(groups))
      
      edges <- data.frame(from = character(), 
                          to = character() , 
                          group = character(), 
                          stringsAsFactors = FALSE)
      
      for (i in iterations){
        temp <- groups[[i]]
        iter <- c(1:nrow(temp))
        for (j in iter){
          from_to <- data.frame(from = temp$TAD[j], to =temp$TAD[1:nrow(temp)] , group = temp$Term[1]) 
          from_to <- from_to[-j,]
          edges <- dplyr::bind_rows(edges, from_to)
        }
      }
      
      edges <- unique(edges)
      g <- graph_from_data_frame(edges, directed = FALSE)
      
      png(filename = paste(string, "/Top 10 ",dir[l]," network graph", ".png", sep = ""), 
          width = 1200, height = 700)
      
      network <- ggraph(g) +
        geom_edge_link(aes(color = group), alpha = 0.5) +     
        geom_node_point(size = 5, shape = 21, stroke = 1,
                        fill = 'white', color = 'black') +
        geom_node_text(aes(label = name), repel = TRUE, size = 4 , fontface = "bold") +
        theme_void()+
        ggtitle(label =paste0("Top 10 ",dir[l]," Network Graph"))+
        theme(
          plot.title = element_text(size=20, hjust = 0.5, vjust = 0.5, face = "bold") ,
          legend.text = element_text(size = 15),
          legend.key.width = unit(1.5,"cm"),
          legend.title = element_text(size = 17)
        )
      
      
      print(network)
      dev.off() 
    }
    
    
    # plot 4
    #group data according to Term 
    if (str_detect(dir[l],"GO")){
      data.visual$ID <- paste0("GO ",data.visual$ID)
      data.plot4 <- data.visual
    }else if (str_detect(dir[l],"KEGG")){
      data.plot4 <- data.visual
      data.plot4$ID <- data.plot4$Term
    }
    
    TAD_grouping <- data.plot4 %>%
      group_by(Term)
    new_grouping <- group_split(TAD_grouping)
    dir.create(paste(string, "/Adjusted P values per ",dir[l], sep = ""), showWarnings = FALSE)
    iterations <- c(1:length(new_grouping))
    
    for (i in iterations) {
      
      temp <- new_grouping[[i]]
      if(l == 3){
        temp$Term <- rep_len("",nrow(temp))
      }
      temp$ID <- str_replace(temp$ID,"/","-")
      png(filename = paste(string, "/Adjusted P values per ",dir[l], "/Adjusted P values of ",temp$ID[1], ".png", sep = ""), 
          #width = 900, height = 500)
          height = (20*nrow(temp)+150), width = 100*ncol(temp))
      p4 <- temp %>%
        ggplot(aes(x = as.factor(reorder(TAD, P.adjust)),y = P.adjust)) +
        geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity", binwidth = 0.5) +
        labs(title = paste0("Adjusted P values of ",temp$ID[1]," in different TADs"), 
             subtitle = temp$Term[1])+
        ylab("Adjusted P value")+
        xlab("TAD number")+
        theme_ipsum() +
        theme(
          plot.title = element_text(size=14, face = "bold"),
          plot.subtitle = element_text(size=10),
          axis.title.x = element_text(size = 13,hjust = 0.5, vjust = 0.5,face = "bold"),
          axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
          
        )+coord_flip() 
      
      print(p4)
      dev.off()
    }
    
    # plot 5
    #group data according to TAD 
    TAD_grouping <- data.visual %>%
      group_by(TAD)
    new_grouping <- group_split(TAD_grouping)
    dir.create(paste(string, "/Adjusted P values per TAD", sep = ""), showWarnings = FALSE)
    iterations <- c(1:length(new_grouping))
    
    for (i in iterations) {
      
      temp <- new_grouping[[i]]
      if (nrow(temp)>30){
        text = "Showing top 30 Terms, for the complete list consult the output csv"
        temp <- temp[1:30,]
      }else{
        text = ""
      }
      png(filename = paste(string, "/Adjusted P values per TAD/","Adjusted P values of ",temp$TAD[1], ".png", sep = ""), 
          #width = 900, height = 500)
          height = (20*nrow(temp)+150), width = 1000)
      p5 <- temp %>%
        ggplot(aes(x = as.factor(reorder(Term, P.adjust)),y = P.adjust)) +
        geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity") +
        labs(title = paste0("Adjusted P values of ",dir[l]," in ",temp$TAD[1]), 
             subtitle = text)+
        ylab("Adjusted P value")+
        xlab(dir[l])+
        theme_ipsum() +
        theme(
          plot.title = element_text(size=14, face = "bold"),
          plot.subtitle = element_text(size=11),
          axis.title.x = element_text(size = 11.5,hjust = 0.5, vjust = 0.5,face = "bold"),
          axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
          
        )+coord_flip() 
      
      print(p5)
      dev.off()
    }
  }
}


pathVisual <- function(full.data, pathview.input , output_folder){
  
  pathview.input <- pathview.input[which(pathview.input$hsa.ids != "NA"),]
  #plot - pathview
  data.S <- full.data %>% 
    dplyr::select(Gene_id, diff)
  data.S <- data.S %>% 
    separate_rows(Gene_id, sep = "\\|") %>%
    as.data.table()
  data.S <- data.S[which((data.S$Gene_id != "NA") &(data.S$Gene_id !="")), ]
  
  new.dir <- paste0(output_folder,"/Pathview")
  dir.create(new.dir, showWarnings = FALSE)
  
  #load again gene ID types because pathview library doesn't load properly
  #after clearing the workspace, so just to be sure
  data(gene.idtype.bods)
  
  loops <- c(1:nrow(pathview.input))
  
  for (l in loops){
    
    path <- pathview.input[l]
    path <- path %>%
      separate_rows(Genes, sep = ";", convert = TRUE) %>%
      unique() %>%
      as.data.table()
    path <- merge(path, data.S,by.x = "Genes", by.y = "Gene_id")
    path <- unique(path)
    path$abs_diff <- abs(path$diff)
    path <- setorder(path,Genes,-abs_diff)
    path <- group_by(path,Genes) %>%
      dplyr::summarise(hsa.ids,Term,Genes,diff = diff[1]) %>%
      as.data.table() %>%
      unique()
    
    path <- path %>% remove_rownames %>% column_to_rownames(var = "Genes")
    
    if (nrow(path)>1){
      
      dir <- paste0(new.dir,"/",path$hsa.ids[1])
      dir.create(dir, showWarnings = FALSE)
      
      current.folder <- paste0(getwd(),"/",path$hsa.ids[1],".pathview.png")
      new.folder <- paste0(dir,"/",path$hsa.ids[1],".pathview.png")
      
      p.input <- as.matrix(path[,3])
      rownames(p.input) <- rownames(path)
     
       p <- pathview(gene.data  = p.input,
                    pathway.id = path$hsa.ids[1],
                    species    = "hsa",
                    gene.idtype = "SYMBOL",
                    same.layer = F,   #two-layer graph (node colors and labels are added + official gene symbols)
                    kegg.dir = dir ,
                    kegg.native = T,
                    na.col = "white",keys.align = "y",multi.state = F) 
      
      
      file.copy(current.folder,new.folder,overwrite = T)
      file.remove(current.folder)
      
    }
  }
}

motifVisual <- function(image_output_folder,motif_output_folder){
  
  report.list <- dget(paste0(motif_output_folder,"/report_motif.txt"))
  data.plot3.4 <- fread(paste0(motif_output_folder, "/data for TFs visualization.csv"))
  
  #plot1+2
  dir.create(paste(image_output_folder, "/Plots per TADs", sep = ""), showWarnings = FALSE)
  iterations <- c(1:length(report.list))
  
  
  data.plot5.6 <- data.table(P.value = numeric(),
                             RawScores = numeric())
  for (i in iterations){
    
    tad.dir = paste0(image_output_folder, "/Plots per TADs/",report.list[[i]]@d$tad[1])
    dir.create(tad.dir,showWarnings = F)
    
    #plot 1
    png(filename = paste(tad.dir, "/PWMEnrich_Image.png", sep = ""),width = 1100, height = 1200)
    
    temp <-report.list[[i]]
    temp@d <- temp@d[,1:6]
    plot(temp[1:10], fontsize=20, id.fontsize=20)
    
    dev.off()
    
    tabl <- report.list[[i]]@d
    
    tableFor5.6 <- tabl %>%
      dplyr::select(p.value,raw.score)
    colnames(tableFor5.6) <- c("P.value", "RawScores")
    
    data.plot5.6 <- rbind(data.plot5.6,tableFor5.6)
    
    #exclude uncharacterized motifs
    tabl <- tabl[str_detect(tabl$id,"UW.Motif.",negate = TRUE),]
  
    if (nrow(tabl)>30){
      text = "Showing top 30 Terms, for the complete list consult the over-represented TFs in each tad.csv"
      tabl <- tabl[1:30,]
    }else{
      text = ""
    }
    
    #plot 2
    png(filename = paste(tad.dir, "/Histogram Image.png", sep = ""), 
        #width = 900, height = 500)
        height = (20*nrow(tabl)+150), width = 1000)
    
    p2 <- tabl %>%
      ggplot(aes(x = as.factor(reorder(target, p.value)),y = p.value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity") +
      labs(title = paste0("P values of TFs in ",tabl$tad[1]), 
           subtitle = text)+
      ylab("P value")+
      xlab("Transcription Factors")+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=14, face = "bold"),
        plot.subtitle = element_text(size=11),
        axis.title.x = element_text(size = 11.5,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 13,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 
    
    print(p2)
    dev.off()
  }
  
  #plot 3+4
  #get PWMs for TFs
  loops <- c(1:length(report.list))
  PWMs <- list()
  for (l in loops){
    temp <- report.list[[l]]@pwms
    PWMs <- c(PWMs,temp)
  }
  PWMs <- unique(PWMs)
  
  loops <-c(1:length(PWMs))
  names.list <- c()
  for(l in loops){
    names.list <- c(names.list,PWMs[[l]]@id)
  }
  
  names(PWMs) <- names.list
  
  dir.create(paste(image_output_folder, "/Transcription Factors", sep = ""), showWarnings = FALSE)
  data.plot3.4 <- data.plot3.4 %>%
    separate_rows(P.value, sep = "\\|")
  data.plot3.4 <- data.plot3.4 %>%
    group_by(TFs,TAD)%>%
    summarise(TFs,TAD, P.value = min(P.value),Motifs) %>%
    unique()
  data.plot3.4 <- data.plot3.4[which(data.plot3.4$TFs != "NA"),]
  data.plot3 <- data.plot3.4 %>%
    group_by(TFs)
  groups.TFs <- group_split(data.plot3) 
  iterations <- c(1:length(groups.TFs))
  
  for (i in iterations) {
    
    temp <- groups.TFs[[i]]
    temp$folder.name <- temp$TFs
    charac <- c(1:round(nchar(temp$folder.name[1])/2))
    for(c in charac){
      temp$folder.name <- str_replace(temp$folder.name,"/","-")
      temp$folder.name <- str_replace(temp$folder.name," / ","-")
      temp$folder.name <- str_replace(temp$folder.name,":","-") 
    }
    string <- paste(image_output_folder, "/Transcription Factors/",temp$folder.name[1], sep = "")
    dir.create(string, showWarnings = F)
    
    #plot3
    temp$P.value <- as.numeric(temp$P.value)
    temp <- setorder(temp, -P.value)
    
    if (nrow(temp)>30){
      text = "Showing top 30 Terms, for the complete list consult the output csv"
      temp <- temp[1:30,]
    }else{
      text = ""
    }
    
    png(filename = paste(string, "/Barplot.png", sep = ""), 
        #width = 900, height = 500)
        height = (20*nrow(temp)+150), width = 1000)
    p3 <- temp %>%
      ggplot(aes(x = as.factor(reorder(TAD,P.value)),y = P.value)) +
      geom_histogram(fill="#66ccff", color="#e9ecef", alpha=0.9, stat = "identity", binwidth = 0.5) +
      labs(title = paste0("P values of ",temp$folder.name[1]," in different TADs"), 
           subtitle = text)+
      ylab("P value")+
      xlab("TAD number")+
      theme_ipsum() +
      theme(
        plot.title = element_text(size=16, face = "bold"),
        plot.subtitle = element_text(size=10),
        axis.title.x = element_text(size = 15,hjust = 0.5, vjust = 0.5,face = "bold"),
        axis.title.y = element_text(size = 15,hjust = 0.5, vjust = 0.5, face = "bold")
        
      )+coord_flip() 
    
    print(p3)
    dev.off()
    
    #plot4
    motifs <- c(str_split(temp$Motifs[1], pattern = "\\|",simplify = T))
    matrices.motifs <- list()
    
    loops <- c(1:length(motifs))
    for (l in loops){
      matrices.motifs[[l]] <- PWMs[[motifs[l]]]@pwm
    }
    names(matrices.motifs) <- motifs
    
      png(filename = paste(string, "/Motifs.png", sep = ""))
      
      p4 <- ggseqlogo(matrices.motifs, method = 'prob')
      
      print(p4)
      dev.off() 
    
  }
  
  #plot5 - P.value and P.adjust Density Plots
  
  png(filename = paste(image_output_folder, "/","Density plot-P values of Motifs.png", sep = ""), 
      width = 500, height = 650)
  
  p5 <- ggplot(data=data.plot5.6, aes(x=P.value)) +
    geom_density(adjust=1.5, fill="#69b3a2", color="#e9ecef", alpha=0.8) +
    ggtitle("P values of Motifs") +
    xlab("P value")+
    ylab("Density")+
    theme_ipsum()+
    theme(
      plot.title = element_text(size=15),
      axis.title.x = element_text(size = 15, vjust = 0.5,hjust = 0.5),
      axis.title.y = element_text(size = 15, vjust = 1.5,hjust = 0.5),
      legend.text = element_text(size=14),
      legend.title = element_text(size=15)
    )
  
  print(p5)
  dev.off()
  
}

