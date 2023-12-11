### load function annotation file
loadanno <- function(term.type = "Pathway", diseaseid, dis.geneid){
  #疾病基因富集分析
  #自行下载kegg 和reactome数据库中 基因-通路关系; GO 基因-本体关系; HPO 基因-表型关系; Drug 基因-药物关系####

  term2gene <- system.file("extdata", paste0(term.type, '_anno.Rdata'), package = "dkgminer")
  load(term2gene)

  ### TODO: 补充不同term对应的不同路径下的Rdata文件   term.gene
  if(term.type == "Pathway" | term.type == "GO"){

    library(clusterProfiler)
    term.enrich.df <-enricher(dis.geneid,
                              TERM2GENE=term.gene,
                              TERM2NAME=NA,
                              qvalueCutoff=1,
                              pvalueCutoff=1,
                              minGSSize=0,
                              maxGSSize = 50000)
    term.enrich<-term.enrich.df@result
    term.enrich<-term.enrich[term.enrich$p.adjust<=0.05,]
  }

  if(term.type == "Pheno" | term.type == "Drug"){

    term2enrich <- system.file("extdata", paste0(term.type, '_enrich.Rdata'), package = "dkgminer")
    load(term2enrich)


    term.enrich <- term.enrich[term.enrich$DiseaseID==diseaseid, ]
  }

  return( list(term.gene = term.gene, term.enrich = term.enrich) )
}





### add common function
addcomfunc <- function(term.type = "Pathway", net.edge, term.gene, term.enrich){
  edge.df.cs.function <- data.frame()
  edge.df.cs2 <- data.frame()

  for(i in 1:dim(net.edge)[1] ){
    # i<-6
    x1<-as.character(net.edge[i,]$start)
    x2<-as.character(net.edge[i,]$end)

    y<-intersect(term.gene[grep(x1,term.gene$gene),]$term, term.gene[grep(x2,term.gene$gene),]$term)


    y2<-y
    if(length(y2)>0){
      enrich.y2<-term.enrich[term.enrich$ID %in% y2,]
      if(term.type == "Pathway" | term.type == "GO"){enrich.y2<-enrich.y2[order(enrich.y2$Count,decreasing = TRUE),]}

      temp.F<-data.frame(term=y2)
      temp.F$pro1<-rep(x1,each=dim(temp.F)[1])
      temp.F$pro2<-rep(x2,each=dim(temp.F)[1])
      temp.F$Type<-rep("only.co-occurrence",each=dim(temp.F)[1])
      if(dim(temp.F[temp.F$term %in% enrich.y2$ID, ])[1]>0){
        temp.F[temp.F$term %in% enrich.y2$ID, ]$Type<-"Enriched"
      }
      temp.F$Enrich.count<-rep(dim(enrich.y2)[1],each=dim(temp.F)[1])
      temp.F$cooccur.count<-rep((length(y)),each=dim(temp.F)[1])
      temp.F$onto.type<-rep(term.type, each=dim(temp.F)[1])

      temp.F<-temp.F[,c(2,3,1,4,5,6,7)]
      edge.df.cs.function<-rbind(edge.df.cs.function,temp.F)
    }

    if(length(y)>0){
      enrich.y<-term.enrich[term.enrich$ID %in% y,]
      if(dim(enrich.y)[1]>=5){
        if(term.type == "Pathway" | term.type == "GO"){y1<-enrich.y[order(enrich.y$Count,decreasing = TRUE),]$ID[1:5]} else {y1<-enrich.y$ID[1:5]}
        y<-y1
      }
      if(dim(enrich.y)[1]>0&dim(enrich.y)[1]<5){
        if(term.type == "Pathway" | term.type == "GO"){y1<-enrich.y[order(enrich.y$Count,decreasing = TRUE),]$ID[1:dim(enrich.y)[1]]} else {y1<-enrich.y$ID[1:dim(enrich.y)[1]]}

        if(length(setdiff(y,y1))>=4){
          y<-c(y1,setdiff(y,y1)[1:(5-length(y1))])
        }else{
          y<-c(y1,setdiff(y,y1))
        }
      }

      if(dim(enrich.y)[1]==0&length(y)>=5){
        y<-y[1:5]
      }
      if(dim(enrich.y)[1]==0&length(y)<5){
        y<-y
      }

      temp<-data.frame(pro1=c(rep(x1,length(y)),y))
      temp$pro2<-c(y,rep(x2,length(y)))
      temp$start<-rep(x1,dim(temp)[1])
      temp$end<-rep(x2,dim(temp)[1])
      temp$method<-rep("common.function",dim(temp)[1])
      edge.df.cs2<-rbind(temp,edge.df.cs2)
    }
  }

  return( list(edge.df = edge.df.cs2, edge.df.function = edge.df.cs.function) )
}





##### 生成网络 边&节点 文件
### candi.gene，用户选择的预测基因
outnet <- function(candi.gene, omics.gene, rw.evidence, cs.evidence, g, term.type = c("Pathway", "GO", "Pheno", "Drug"), diseaseid, dis.cutoff = 0.3, ge.method = "deepwalk"){
  ### process # disease genes
  disease.gene.file <- system.file("extdata", "Disease_associated_Gene.csv", package = "dkgminer")
  disease.gene <- fread(disease.gene.file)

  dis <- disease.gene[disease.gene$DiseaseId==diseaseid & disease.gene$score>dis.cutoff, ]

  dis.gene<-dis%>%
    dplyr::group_by(DiseaseId,GeneId)%>%
    dplyr::summarise(DiseaseId=unique(DiseaseId),
                     GeneId=unique(GeneId),
                     score=max(score,na.rm = TRUE))

  dis.geneid <- unique(dis.gene$GeneId)
  ###



  #1.对于随机行走方法得到的证据对，采用最短路径的方式来呈现两蛋白之间关系
  evidence.df<-rw.evidence[rw.evidence$pre.genes == candi.gene,]
  edge.df.rw <- data.frame()

  if(dim(evidence.df)[1]>0){
    for(i in 1:dim(evidence.df)[1]){
      # i<-3
      x<-evidence.df$pre.genes[i]
      y<-evidence.df$known.dis.genes[i]
      results1 <-shortest_paths(g, from = V(g)[name == x], to =  V(g)[name==y], output = 'epath')$epath
      if(length(results1[[1]])>0){
        results1 <- sapply(results1, as_ids)

        temp <- as.data.frame(matrix(results1))
        temp$start<-rep(x,dim(temp)[1])
        temp$end<-rep(y,dim(temp)[1])
        edge.df.rw <- rbind(temp,edge.df.rw)
      }

    }
    if(dim(edge.df.rw)[1]>0){
      edge.df.rw<-separate(edge.df.rw,V1,c("pro1","pro2"),"\\|")
      # edge.df.rw$pro1.type<-rep("overlap",dim(edge.df.rw.only)[1])
      # edge.df.rw$pro2.type<-rep("known.dis.rw",dim(edge.df.rw.only)[1])
      edge.df.rw$method<-rep("rw",dim(edge.df.rw)[1])
    }
  }



  # 2.对于向量嵌入方法得到的证据对，采用直接关联方式来呈现两蛋白之间关系
  evidence.df2<-cs.evidence[cs.evidence$pre.genes == candi.gene,]

  if(dim(evidence.df2)[1]>0){
    if(TRUE){
      edge.df.cs<-data.frame(pro1=evidence.df2$pre.genes,
                             pro2=evidence.df2$known.dis.genes)
      edge.df.cs$start<-edge.df.cs$pro1
      edge.df.cs$end<-edge.df.cs$pro2
      # edge.df.cs$method<-rep("deepwalk",dim(edge.df.cs)[1])
      edge.df.cs$method<-rep(ge.method, dim(edge.df.cs)[1])
    }
  }



  # 3. 为蛋白关系对，寻找功能关联
  if(dim(edge.df.rw)[1] > 0 & dim(edge.df.cs)[1] > 0){
    edge.df.F<-unique(rbind(edge.df.rw[,c("start","end")],edge.df.cs[,c("start","end")]))
  }
  if(dim(edge.df.rw)[1] > 0 & dim(edge.df.cs)[1] == 0){
    edge.df.F<-unique(rbind(edge.df.rw[,c("start","end")]))
  }
  if(dim(edge.df.rw)[1] == 0 & dim(edge.df.cs)[1] > 0){
    edge.df.F<-unique(rbind(edge.df.cs[,c("start","end")]))
  }

  # edge.df.cs2<-data.frame()
  # edge.df.cs.function<-data.frame()

  ### 根据用户选择的本体， 增加 共有功能的边
  if(length(term.type) > 1){
    edgeresult.ls <- list()
    for(i in 1:length(term.type)){
      terminfo <- loadanno(term.type[i], diseaseid, dis.geneid)
      edgeresult.ls[[i]] <- addcomfunc(term.type[i], net.edge = edge.df.F, terminfo$term.gene, terminfo$term.enrich)
    }
    for(i in 1:(length(term.type)-1)){
      y <- i + 1
      edgeresult.ls[[y]]$edge.df <- rbind(edgeresult.ls[[i]]$edge.df, edgeresult.ls[[y]]$edge.df)
      edgeresult.ls[[y]]$edge.df.function <- rbind(edgeresult.ls[[i]]$edge.df.function, edgeresult.ls[[y]]$edge.df.function)
    }

    edgeresult <- edgeresult.ls[[length(term.type)]]
  } else {
    terminfo <- loadanno(term.type, diseaseid, dis.geneid)
    edgeresult <- addcomfunc(term.Type, net.edge = edge.df.F, terminfo$term.gene, terminfo$term.enrich)
  }





  edge.df.cs.function <- edgeresult$edge.df.function
  edge.df.cs2 <- edgeresult$edge.df

  edge.df.cs.function$pro.pair<-paste(edge.df.cs.function$pro1,edge.df.cs.function$pro2,sep=",")
  edge.df.cs.function[edge.df.cs.function$pro1>edge.df.cs.function$pro2,]$pro.pair<-paste(edge.df.cs.function[edge.df.cs.function$pro1>edge.df.cs.function$pro2,]$pro2,edge.df.cs.function[edge.df.cs.function$pro1>edge.df.cs.function$pro2,]$pro1,sep=",")

  edge.df<-rbind(edge.df.rw,edge.df.cs)
  edge.df<-rbind(edge.df,edge.df.cs2)

  edge.df$pro.pair<-paste(edge.df$pro1,edge.df$pro2,sep=",")
  if(dim(edge.df[edge.df$pro1>edge.df$pro2,])[1]>0){
    edge.df[edge.df$pro1>edge.df$pro2,]$pro.pair<-paste(edge.df[edge.df$pro1>edge.df$pro2,]$pro2,edge.df[edge.df$pro1>edge.df$pro2,]$pro1,sep=",")
  }
  edge.df<-merge(edge.df,unique(edge.df.cs.function[,-c(1:4)]), all.x=TRUE)

  edge.df<-as.data.frame(separate(edge.df,pro.pair,c("pro1","pro2"),remove = FALSE,sep=","))
  edge.df<-edge.df[, c(2:ncol(edge.df), 1)]



  # 4.扩充边-增加网络中节点的直接关联
  if(TRUE){
    node.df<-data.frame(pro=unique(c(edge.df$pro1, edge.df$pro2)))
    node.df$pro.type<-rep("linker",dim(node.df)[1])
    # node.df[node.df$pro %in% dis.geneid  ,]$pro.type<-"known.dis"


    ### 确定各种蛋白的类型（rw发现的、cs发现的、组学已知的）
    ### TODO: 添加 omics.rank.df1 的来源
    for(i in 1:length(term.type)){
      terminfo <- loadanno(term.type[i], diseaseid, dis.geneid)

      if( nrow(node.df[node.df$pro %in% unique(terminfo$term.gene$term), ]) > 0 ){
        node.df[node.df$pro %in% unique(terminfo$term.gene$term), ]$pro.type <- term.type[i]
      }
      if( nrow(node.df[node.df$pro %in% unique(terminfo$term.enrich$ID), ]) > 0){
        node.df[node.df$pro %in% unique(terminfo$term.enrich$ID),]$pro.type <- paste0("enrich.", term.type[i])
      }
    }


    node.df[node.df$pro %in% candi.gene,]$pro.type<-"candi.gene"
    if(dim(evidence.df)[1]>0){
      node.df[node.df$pro %in% unique(evidence.df$known.dis.genes) ,]$pro.type<-"rw.known.dis"
    }
    if(dim(evidence.df2)[1]>0){
      node.df[node.df$pro %in% unique(evidence.df2$known.dis.genes) ,]$pro.type<-paste0(ge.method, ".known.dis")
    }
    if(dim(evidence.df2)[1]>0&dim(evidence.df)[1]>0){
      node.df[node.df$pro %in% intersect(unique(evidence.df$known.dis.genes),unique(evidence.df2$known.dis.genes)) ,]$pro.type<-paste0("rw.", ge.method, ".known.dis")
    }


    # if(dim(node.df[node.df$pro %in% unique(omics.rank.df1$Genes) ,])[1]>0){
    #   node.df[node.df$pro %in% unique(omics.rank.df1$Genes) ,]$pro.type<-"omics.known.dis"
    # }


    ### Gene ID to Gene Symbol
    node.df$pro.Symbol <- node.df$pro
    node.df[grep("^[0-9]",node.df$pro),]$pro.Symbol <- mapIds(org.Hs.eg.db,
                                                              keys= node.df[grep("^[0-9]",node.df$pro),]$pro,
                                                              column="SYMBOL",
                                                              keytype="ENTREZID",
                                                              multiVals="first")

    #   ### ontology names   TODO: 增加各种本体的信息表格
    #   if(dim(node.df[grep("^path",node.df$pro),])[1]){
    #     node.df[grep("^path",node.df$pro),]$pro.Symbol<-keggid.name[keggid.name$pathway %in% node.df[grep("^path",node.df$pro),]$pro,]$path.name
    #   }

    #   if(dim(node.df[grep("^R",node.df$pro),])[1]){
    #     node.df[grep("^R",node.df$pro),]$pro.Symbol<-unique(reactom.df[,c(2,3)])[unique(reactom.df[,c(2,3)])$pathway %in% node.df[grep("^R",node.df$pro),]$pro,]$path.name
    #   }
    #   if(dim(node.df[grep("^GO",node.df$pro),])[1]){
    #     node.df[grep("^GO",node.df$pro),]$pro.Symbol<-node.df[grep("^GO",node.df$pro),]$pro
    #   }

    #   if(dim(node.df[grep("^DB",node.df$pro),])[1]){
    #     node.df[grep("^DB",node.df$pro),]$pro.Symbol<- drug.info[drug.info$DrugBank_ID %in% node.df[grep("^DB",node.df$pro),]$pro,]$name
    #   }

    #   if(dim(node.df[grep("^HP",node.df$pro),])[1]){
    #     node.df[grep("^HP",node.df$pro),]$pro.Symbol<- pheno.info[pheno.info$ID %in% node.df[grep("^HP",node.df$pro),]$pro,]$Term
    #   }

    #   node.df<-merge(node.df,omics.df,by.x="pro",by.y="Genes",all.x=TRUE)
  }


  ### TODO: 对于网络 边文件 进行筛选，确定最终网络图上呈现的结果
  edge.df1<-unique(edge.df[,-c(1:4)])
  # edge.df2<-edge.df1[edge.df1$method!="deepwalk", ]
  edge.df2<-edge.df1[edge.df1$method!=ge.method, ]
  edge.df3<-edge.df2 %>%
    dplyr::group_by(pro.pair)%>%
    dplyr::summarise(
      method=paste(unique(method),collapse = ","),
      Enrich.count=sum(Enrich.count),
      cooccur.count=sum(cooccur.count),
      onto.type=paste(unique(onto.type),collapse = ",")
    )

  # edge.df22<-edge.df1[edge.df1$method=="deepwalk", ]
  edge.df22<-edge.df1[edge.df1$method==ge.method, ]
  edge.df33<-edge.df22 %>%
    dplyr::group_by(pro.pair)%>%
    dplyr::summarise(
      method=paste(unique(method),collapse = ","),
      Enrich.count=sum(Enrich.count),
      cooccur.count=sum(cooccur.count),
      onto.type=paste(unique(onto.type),collapse = ",")
    )
  edge.df3<-rbind(edge.df3,edge.df33)
  edge.df3<-as.data.frame(separate(edge.df3,pro.pair,c("pro1","pro2"),remove = FALSE,sep=","))
  edge.df4<- edge.df3[,c(2,3,4,5,6,7,1)]
  edge.df4[edge.df4$onto.type=="NA",]$onto.type<-""








  ### node informations
  # 功能注释相关 节点信息
  node.anno1<-edge.df.cs.function[,c(-1,-5,-6)]
  names(node.anno1)[1]<-"pro"
  node.anno2<-edge.df.cs.function[,c(-2,-5,-6)]
  names(node.anno2)[1]<-"pro"

  node.anno<-unique(rbind(node.anno1,node.anno2))
  node.anno<-merge(node.anno,node.df[,c(1,3)],by.x="pro",by.y="pro",all.x=TRUE)

  node.anno.G<-node.anno %>%
    dplyr::group_by(pro, onto.type, Type)%>%
    dplyr::summarise(
      term=paste(unique(term), collapse = ";"),
      term.n=n()
    )

  # 随机游走（PPI）相关 节点信息
  node.anno.rw1<-edge.df[edge.df$method=="rw",][,c(1,2,9)]
  node.anno.rw2<-edge.df[edge.df$method=="rw",][,c(2,1,9)]
  names(node.anno.rw2)<-names(node.anno.rw1)
  node.anno.rw<-unique(rbind(node.anno.rw1,node.anno.rw2))

  node.anno.rw$pro1
  node.anno.rw<-merge(node.anno.rw,node.df[,c(1,3)],by.x="pro1",by.y="pro",all.x=TRUE)
  names(node.anno.rw)[names(node.anno.rw)=="pro.Symbol"]<-"pro1.Symbol"
  node.anno.rw<-merge(node.anno.rw,node.df[,c(1,3)],by.x="pro2",by.y="pro",all.x=TRUE)
  names(node.anno.rw)[names(node.anno.rw)=="pro.Symbol"]<-"pro2.Symbol"

  node.anno.rw.G<-node.anno.rw %>%
    dplyr::group_by(pro1)%>%
    dplyr::summarise(
      term=paste(unique(pro2.Symbol), collapse = ";"),
      term.n=n()
    )

  node.anno.rw.G$onto.type<-rep("gene.actors",dim(node.anno.rw.G)[1])
  node.anno.rw.G$Type<-rep("",dim(node.anno.rw.G)[1])
  # node.anno.rw.G<-merge(node.anno.rw.G,node.df[,c(1,2)],by.x="pro1",by.y="pro",all.x=TRUE)
  node.anno.rw.G<-node.anno.rw.G[,c(1,4,5,2,3)]
  names(node.anno.rw.G)<-names(node.anno.G)


  # 节点信息合并
  node.anno.df<-rbind(node.anno.rw.G, node.anno.G)
  node.anno.df$Gene.name<-mapIds(org.Hs.eg.db,
                                 keys=node.anno.df$pro,
                                 keytype="ENTREZID",
                                 column="GENENAME",
                                 multiVals="first")

  node.anno.df$is.known.dis.gene<-rep("No",dim(node.anno.df)[1])
  node.anno.df[node.anno.df$pro %in% dis.geneid,]$is.known.dis.gene<-"Yes"

  node.anno.df$is.omics.gene<-rep("No",dim(node.anno.df)[1])
  node.anno.df[node.anno.df$pro %in% omics.gene,]$is.omics.gene<-"Yes"





  ### term informations
  node.anno.term1<-(edge.df.cs.function[,c(3,1,4)])
  node.anno.term2<-(edge.df.cs.function[,c(3,2,4)])
  names(node.anno.term2)<-names(node.anno.term1)
  node.anno.df2<-unique(rbind(node.anno.term2,node.anno.term1))
  node.anno.df2$pro.sym<-mapIds(org.Hs.eg.db,
                                keys=node.anno.df2$pro1,
                                keytype="ENTREZID",
                                column="SYMBOL",
                                multiVals="first")

  node.anno.df2.G<-node.anno.df2 %>%
    dplyr::group_by(term)%>%
    dplyr::summarise(
      genes=paste((pro.sym), collapse = ";"),
      term.n=n(),
      term.name=unique(term),
      Type=unique(Type)
    )

  node.anno.df2.G$term.type<-"unsure"

  if(dim(node.anno.df2.G[grep("^DB",node.anno.df2.G$term),])[1]>0){
    node.anno.df2.G[grep("^DB",node.anno.df2.G$term),]$term.type<-"Drug"
  }
  if(dim(node.anno.df2.G[grep("^GO",node.anno.df2.G$term),])[1]>0){
    node.anno.df2.G[grep("^GO",node.anno.df2.G$term),]$term.type<-"GO"
  }
  if(dim(node.anno.df2.G[grep("^HP",node.anno.df2.G$term),])[1]>0){
    node.anno.df2.G[grep("^HP",node.anno.df2.G$term),]$term.type<-"Pheno"
  }
  if(dim(node.anno.df2.G[grep("^path",node.anno.df2.G$term),])[1]>0){
    node.anno.df2.G[grep("^path",node.anno.df2.G$term),]$term.type<-"Pathway"
  }
  if(dim(node.anno.df2.G[grep("^R",node.anno.df2.G$term),])[1]>0){
    node.anno.df2.G[grep("^R",node.anno.df2.G$term),]$term.type<-"Pathway"
  }





  ### edge informations
  edge.anno.df <- edge.df.cs.function


  return(list(nodes = node.df, edges = edge.df4, funcls = list(node.anno = node.anno.df, term.anno = node.anno.df2.G, edge.anno = edge.anno.df) ))

}





plotvisnet <- function(nodes, edges, annoinfo){
  ### 1108 添加注释信息
  # node anno info
  annoinfo_node <- annoinfo$node.anno

  scrollbox <- '<div style="width:400px;height:300px;line-height:2em;overflow:scroll;padding:1px;">'

  titletext <- c()
  for(i in 1:length(unique(annoinfo_node$pro))){
    nn <- which(annoinfo_node$pro == unique(annoinfo_node$pro)[i])
    textls <- lapply(nn, function(x) paste0("<br>List of ", annoinfo_node$onto.type[x], "[", annoinfo_node$Type[x], "]", "(", annoinfo_node$term.n[x], ") :<br>",
                                            "&nbsp;&nbsp;", annoinfo_node$term[x]) )
    textls <- gsub(', ', '<br>', gsub('[]', '', toString(unlist(textls)), fixed = T))
    textls <- paste0("Gene description: ", unique(annoinfo_node$Gene.name[nn]),
                     "<br>Query gene: ", unique(annoinfo_node$is.omics.gene[nn]),
                     "<br>Known disease associated gene: ", unique(annoinfo_node$is.known.dis.gene[nn]),
                     textls)

    titletext[i] <- paste0(scrollbox, textls, '</div>')
  }

  gsym <- mapIds(org.Hs.eg.db,
                 keys= unique(annoinfo_node$pro),
                 column="SYMBOL",
                 keytype="ENTREZID",
                 multiVals="first")

  annoinfo_node <- cbind(pro=gsym, title=titletext)


  # term anno info
  annoinfo_term <- annoinfo$term.anno
  annoinfo_term <- na.omit(merge(nodes, annoinfo_term, by.x='pro', by.y='term', all.x=T))

  annoinfo_term$title <- ''
  annoinfo_term$title <- paste0(annoinfo_term$term.type, ": ", annoinfo_term$term.name, "<br>",
                                "Disease associated: No", "<br>",
                                "List of associated genes(", annoinfo_term$term.n, "):<br>",
                                "&nbsp;&nbsp;", annoinfo_term$genes)

  annoinfo_term <- annoinfo_term[, c("pro", "title")]


  # merge info
  annoinfotab <- rbind(annoinfo_node, annoinfo_term)




  ### node & edge plot
  selgene <- nodes$pro.Symbol[which(nodes$pro.type == 'candi.gene')]

  ### 11.18 temp
  nn <- union(union(grep('known.dis', nodes$pro.type), grep('linker', nodes$pro.type)), grep('candi.gene', nodes$pro.type))
  nodes <- nodes[order(nodes[nn, 'pro.type']), ]
  edges <- edges[-which(edges$method == 'common.function'), ]


  mapvs <- cbind(id=1:nrow(nodes), pro=nodes$pro)

  nodes <- cbind(id=1:nrow(nodes), nodes[, c(3, 2)])
  colnames(nodes) <- c('id', 'label', 'pro.type')

  ### 合并 添加节点注释信息
  # nodes <- cbind(nodes, title=nodes$label)
  nodes <- merge(nodes, annoinfotab, by.x='label', by.y='pro', all = T)

  nn <- which(nchar(nodes$label) > 10)
  nodes$label[nn] <- paste0(substr(nodes$label[nn], start = 1, stop = 10), '...')

  edges <- edges[, 1:5]
  colnames(edges)[1:3] <- c('from', 'to', 'title')
  edges$from <- mapvalues(edges$from, from = mapvs[, 2], to = mapvs[, 1], warn_missing = F)
  edges$to <- mapvalues(edges$to, from = mapvs[, 2], to = mapvs[, 1], warn_missing = F)



  ### nodes info
  nodes$color.background <- 'grey'
  nodes$color.background[grep('linker', nodes$pro.type)] <- 'pink'
  nodes$color.background[grep('enrich.', nodes$pro.type)] <- '#9999FF'
  nodes$color.background[which(nodes$pro.type == 'candi.gene')] <- 'red'   #'orange'
  nodes$color.background[grep('known.dis', nodes$pro.type)] <- '#FFFFFF'

  # https://github.com/datastorm-open/visNetwork/issues/51
  nodes$color.border <- ''
  nodes$color.border[grep('known.dis', nodes$pro.type)] <- 'red'

  nodes$shape <- 'dot'
  nodes$shape[grep('known.dis', nodes$pro.type)] <- 'box'
  nodes$shape[grep('candi.gene', nodes$pro.type)] <- 'ellipse'
  nodes$shape[grep('linker', nodes$pro.type)] <- 'box'

  # https://rdrr.io/cran/shp2graph/man/Degree.list.html
  # nodes$value <- shp2graph::Degree.list(nodes[, 1:2], edges[, 1:2], Directed = F)[[1]]

  ### 通过比例缩放设定每个节点的大小： https://github.com/datastorm-open/visNetwork/issues/49
  nodes$value <- 1
  nodes$value[grep('known.dis', nodes$pro.type)] <- 1.2
  nodes$value[grep('candi.gene', nodes$pro.type)] <- 1.5
  nodes$value[grep('linker', nodes$pro.type)] <- 1

  # nodes$font.size <- 8
  nodes$shadow <- FALSE



  ### edges info
  edges$dashes <- TRUE
  edges$dashes[which(edges$title == 'rw')] <- FALSE

  # edges$color <- 'grey'
  # edges$color[grep('rw', edges$title)] <- 'blue'
  # edges$color[grep('deepwalk', edges$title)] <- 'blue'
  edges$color <- 'blue'
  edges$color[grep('common.function', edges$title)] <- 'grey'

  edges$smooth <- TRUE
  edges$smooth[grep('common.function', edges$title)] <- FALSE

  ### TODO 0727
  edges$label <- ''

  # simN <- grep('deepwalk', edges$title)
  # simN <- which(edges$title != 'common.function' & edges$title != 'rw')
  simN <- which(edges$title != 'common.function')
  # repl <- NULL
  # for(x in 1:length(simN)) {repl <- c(repl, sample(1:3, 1, replace = F))}

  edges$label[simN] <- edges$cooccur.count[simN]


  edges$font.color <- ''
  edges$font.color[simN] <- 'red'





  mt <- list(text = paste0('Network view  for selected  candidate gene : ', selgene),
             style = 'font-family:Georgia, "Times New Roman", Times, serif;color:#000000;font-size:18px;text-align:center;')
  st <- list(text = "Clicking each node/edge to view a pop-up window showing the detailed information",
             style = "font-family:Comic Sans MS;color:#ff0000;font-size:12px;text-align:center;")


  # ### TODO 0803 click information
  # clickfile <- paste0('D:/Project/R_project/RpackTest/netdat/dkgdata/node', s, '_anno.csv')
  # clickinfo <- read.table(clickfile, sep = ',', header = T)

  # clickinfo <- 'Gene description: Re1 silencing transcription factor\n
  #     Query gene: Yes\n
  #     Known disease associated gene:  No\n
  #     List of interactors (1):\n
  #         AKT2, AKT serine/threonine kinase 2 (TF-TG)\n
  #     List of pathways (2)\n
  #         Transcriptional misregulation in cancer\n
  #         p53 signaling pathway\n
  #     List of phenotypes (2)\n
  #         a\n
  #         b\n
  #     List of Drugs (3)\n
  #         c\n
  #         d
  #     '
  # # 0803 test
  # nodes$title <- clickinfo

  ### 1008 passing custom nodes and/or edges
  lnodes <- data.frame(label = c("know disease genes", "linker", "candidate disease gene", "co-occurring functions", "enriched functions"),
                       shape = c( "dot", "dot", "dot", "triangle", "triangle"),
                       color.background = c("#FFFFFF", "pink", "red", "grey", "#9999FF"),
                       color.border = c("red", "", "", "", ""),
                       title = "Informations")

  ledges <- data.frame(label = c("evidence", "commom function"),
                       color = c("blue", "blue"),
                       dashes = c(FALSE, TRUE),
                       # arrows = c("to", "from"),
                       font.align = "top")



  ### 11.18 temp
  nodes <- na.omit(nodes)


  visNetwork(nodes, edges, width = "100%", main = mt, submain = st) %>%
    ### 通过比例缩放设定每个节点的大小： https://github.com/datastorm-open/visNetwork/issues/49
    ### you can set the size of label (and so of shape) by setting the scaling option
    visNodes( scaling = list(label = list(enabled = T, min = 8, max = 15)) ) %>%
    # visEdges(dashes = edges$dashes, smooth = edges$smooth) %>%
    # 参考原代码 https://rdrr.io/cran/visNetwork/src/R/visExport.R
    visExport(type = "png", label = 'PNG', style= "color: #ffffff;background: #3C90BC;border:none;") %>%
    visOptions(manipulation = FALSE, highlightNearest = TRUE, nodesIdSelection = list(enabled = FALSE)) %>%

    ### 1008 add legend
    visLegend(addNodes = lnodes, addEdges = ledges, position = "right", useGroups = FALSE, ncol = 1) %>%
    ### 1110 增加拖拽的交互效果
    visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
    ### 1110 单独拖拽一个节点其他不动
    visPhysics(maxVelocity = 0.1)

    # https://visjs.github.io/vis-network/examples/network/other/popups.html
    # visEvents(type = "once", select = "function() {alert('first selection');}") %>%
    # visEvents(select = "function() {alert( 'Gene description: Re1 silencing transcription factor\\nQuery gene: YesKnown disease associated gene:  No\\nList of interactors (1):\\n    AKT2, AKT serine/threonine kinase 2 (TF-TG)\\nList of pathways (2):\\n    Transcriptional misregulation in cancer\\n    p53 signaling pathway\\nList of phenotypes (2):\\n    a\\n    b\\nList of Drugs (3):\\n    c\\n    d' );}")

}
