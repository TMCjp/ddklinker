# https://rpubs.com/liam/RpackagesDevelopment

csfunc <- function(ge.method = "deepwalk", cs.input, dis.train, N.top = 5){
  nn <- na.omit(match(dis.train, rownames(cs.input)))
  CStab <- cs.input[nn, -nn]
  CStab$V1<-row.names(CStab)
  CStab.raw<-gather(CStab, NodeNames, Score, -V1)
  CStab<-CStab[,-c(grep("V1",names(CStab)))]

  if(TRUE){

    CS.df.evidence<-CStab.raw%>%
      dplyr::group_by(NodeNames)%>%
      dplyr::top_n(n=N.top,wt=Score)
    names(CS.df.evidence)<-c("known.dis.genes","pre.genes","Score")
    #write.csv(CS.df.evidence,paste("result",gsub("\\:",".",diseaseid),paste("deepwalk.CS.evidence",".dis.cutoff.",dis.cutoff,".N.top.",N.top,".csv",sep = ""),sep = "/"),row.names = FALSE)

    # str(CStab)
    CSms <- apply(CStab, 2, function(x) mean(x[order(x,decreasing=TRUE)[1:N.top]],na.rm=TRUE))

    ### 生成结果矩阵
    CS.out <- as.data.frame(matrix(NA, nrow = length(CSms), ncol = 2))
    CS.out[, 1] <- names(CSms)
    CS.out[, 2] <- CSms
    names(CS.out)<-c("Genes","CS.Score")

    CS.out<-CS.out%>%
      dplyr::group_by(Genes)%>%
      dplyr::summarise(CS.Score=max(CS.Score,na.rm = TRUE))

    if(dim(CS.out)[1]>0){
      # CS.out<-CS.out[CS.out$Score>1E-5,]
      CS.out<-CS.out[order(CS.out$CS.Score,decreasing = TRUE),]
      CS.out$CS.rank<-c(1:dim(CS.out)[1])
      rank.test4<-CS.out
    }

  }

  rank.test.CS<-rank.test4
  names(rank.test.CS)[2:3]<-paste(ge.method, names(rank.test.CS)[2:3], sep = ".")


  return(list(cs.rank = rank.test.CS, cs.evidence = CS.df.evidence))
}






predictFC <- function(diseaseid, dis.cutoff = 0.3, bgNet, ge.method = "deepwalk", actor.cutoff = 0, N.top = 5, top.cutoff = 20){
  library(data.table)
  library(tidyr)
  library(igraph)
  library(graph)
  library(RandomWalkRestartMH)


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



  ######################
  #1. 向量嵌入-准备####
  ######################
  ### 读入余弦值打分结果矩阵
  gefile <- system.file("extdata", paste0(bgNet, '_ge_', ge.method, "_cosine.csv"), package = "dkgminer")
  CS.input <- as.data.frame(fread(gefile, sep = ',', header = T))
  row.names(CS.input)<-gsub("\n","",as.character(CS.input[,1]))
  CS.input<-CS.input[,-1]


  ######################
  #2. 随机行走方法-准备####
  ######################
  # Random.time<-1000
  ppinet <- system.file("extdata", paste0(bgNet, '_network.Rdata'), package = "dkgminer")
  load(ppinet)

  ###
  ### 增加判断，igraph 3.5以上的版本，igraph对象略有不同，需要单独处理
  vnum <- as.numeric(substr(igraph.version(), 3, 5))
  if(vnum >= 3.5){
    PPI_MultiplexObject <- create.multiplex(list(PPI_Network))
  } else {
    PPI_MultiplexObject <- create.multiplex(PPI_Network)
  }

  AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
  AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
  ###

  subnet <- as_graphnel(PPI_Network)
  net.edges<-adj(subnet, nodes(subnet))
  net.actor<-data.frame(genes=names(net.edges))
  net.actor$N.actors<-apply(net.actor, 1, function(x) length(unlist(net.edges[x])))


  SeedGene<-intersect(dis.geneid, V(PPI_Network))
  SeedGene.actor<-net.actor[net.actor$genes %in% SeedGene,]
  SeedGene<-SeedGene.actor[SeedGene.actor$N.actors>actor.cutoff, ]$genes
  otherGene<-setdiff(V(PPI_Network), dis.geneid)



  if(length(SeedGene)>0){

    # print(N.top)
    dis.train<-SeedGene



    ### out tab
    rank.out<-data.frame()


    #################
    #3.1 运行随机行走####
    #按基因set作为输入，得到基因评分的平均值
    if(TRUE){

      aa.df.raw<-data.frame()
      for(dis.train.gene in dis.train){
        ## We launch the algorithm with the default parameters (See details on manual)
        RWR_PPI_Results <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,
                                                         PPI_MultiplexObject,dis.train.gene,r=0.7)
        # We display the results
        aa.temp<-RWR_PPI_Results$RWRM_Results
        aa.temp$dis.train.gene<-dis.train.gene
        aa.df.raw<-rbind(aa.temp,aa.df.raw)
      }

      aa.df.raw<-aa.df.raw[!(aa.df.raw$NodeNames %in% dis.train), ]
      aa.df<-as.data.frame(spread(aa.df.raw,NodeNames,Score))

      aa.df.evidence<-aa.df.raw%>%
        dplyr::group_by(NodeNames)%>%
        dplyr::top_n(n=N.top,wt=Score)
      names(aa.df.evidence)<-c("pre.genes","Score","known.dis.genes")

      rw.evidence <- aa.df.evidence
      # write.csv(aa.df.evidence,paste("result",gsub("\\:",".",diseaseid),paste("RW.evidence",".dis.cutoff.",dis.cutoff,".N.top.",N.top,".csv",sep = ""),sep = "/"),row.names = FALSE)


      aa.df.new <- apply(aa.df[,-1], 2, function(x) mean(x[order(x,decreasing=TRUE)[1:N.top]]) )
      aa.df.new <-as.data.frame(aa.df.new)
      names(aa.df.new)<-"RW.Score"
      aa.df.new$Genes<-row.names(aa.df.new)

      aa<-aa.df.new[,c(2,1)]

      if(dim(aa)[1]>0){
        # aa<-aa[aa$Score>1E-5,]
        aa<-aa[order(aa$RW.Score,decreasing = TRUE),]
        aa$RW.rank<-c(1:dim(aa)[1])
        rank.test2<-aa
        # names(rank.test2)[c(2,3)]<-paste( names(rank.test2)[c(2,3)],N.top,sep=".")
      }

      rank.test.RW<-rank.test2
    }



    #################
    #3.2 运行向量嵌入####
    if(TRUE){

      rank.test.CS<-data.frame()

      ### CS.input需要在上面开头定义，根据嵌入方法选择不同的相似性结果文件
      test.CS <- csfunc(ge.method = ge.method, cs.input = CS.input, dis.train = dis.train, N.top = N.top)

      rank.test.CS <- test.CS$cs.rank
      cs.evidence <- test.CS$cs.evidence

    }


    rank.out<-merge(rank.test.CS,rank.test.RW,all=TRUE)
    rank.out$dis<-rep(diseaseid, each=dim(rank.out)[1])
    rank.out$Ntop<-rep(N.top,each=dim(rank.out)[1])
    rank.out$dis.cutoff<-rep(dis.cutoff,each=dim(rank.out)[1])
    #   write.csv(rank.out,paste("result",gsub("\\:",".",diseaseid),paste("rank.df.dis.cutoff.",dis.cutoff,".N.top",N.top,".csv",sep = ""),sep = "/"),row.names = FALSE)



    ### 排序及打分 标准化
    nn <- match(paste0(ge.method, '.CS.rank'), colnames(rank.out))

    rank.out$Rank<-rank.out[, nn]
    rank.out[rank.out[, nn] > rank.out$RW.rank & !is.na(rank.out$RW.rank), ]$Rank <-
      rank.out[rank.out[, nn] > rank.out$RW.rank & !is.na(rank.out$RW.rank), ]$RW.rank
    rank.out$Rank.M <- nrow(rank.out)-rank.out$Rank+1

    rank.out<-rank.out[order(rank.out$Rank.M,decreasing = TRUE),]
    rank.out.1<-rank.out[1:100,]
    rank.out.1$Score<-0.9+(0.999-0.9)/(max(rank.out.1$Rank.M)-min(rank.out.1$Rank.M))*(rank.out.1$Rank.M-min(rank.out.1$Rank.M))

    rank.out.2<-rank.out[101:200,]
    rank.out.2$Score<-0.8+(0.9-0.8)/(max(rank.out.2$Rank.M)-min(rank.out.2$Rank.M))*(rank.out.2$Rank.M-min(rank.out.2$Rank.M))

    rank.out.3<-rank.out[201:500,]
    rank.out.3$Score<-0.6+(0.8-0.6)/(max(rank.out.3$Rank.M)-min(rank.out.3$Rank.M))*(rank.out.3$Rank.M-min(rank.out.3$Rank.M))

    rank.out.4<-rank.out[500:dim(rank.out)[1],]
    rank.out.4$Score<-0.01+(0.6-0.01)/(max(rank.out.4$Rank.M)-min(rank.out.4$Rank.M))*(rank.out.4$Rank.M-min(rank.out.4$Rank.M))

    rank.out<-rbind(rank.out.1,rank.out.2,rank.out.3,rank.out.4)
    rank.out$sup.algo<-rep("low.score",dim(rank.out)[1])
    if(dim(rank.out[is.na(rank.out$Score),])[1]>0){
      rank.out[is.na(rank.out$Score),]$sup.algo<-"no.score"
    }

    if(length(which( (rank.out[, nn] <= top.cutoff) == "TRUE" )) > 0 ){ rank.out[rank.out[, nn] <= top.cutoff, ]$sup.algo <- ge.method }
    if(length(which( (rank.out$RW.rank <= top.cutoff & !is.na(rank.out$RW.rank)) == "TRUE" )) > 0 ){ rank.out[rank.out$RW.rank <= top.cutoff & !is.na(rank.out$RW.rank), ]$sup.algo <- "rw" }
    if(length(which( (rank.out[, nn] <= top.cutoff & rank.out$RW.rank <= top.cutoff) == "TRUE" )) > 0 ){ rank.out[rank.out[, nn] <= top.cutoff & rank.out$RW.rank <= top.cutoff, ]$sup.algo <- paste0("rw.", ge.method) }

    rank.out <- rank.out[order(rank.out$Score,decreasing = TRUE), ]
  }


  return(list(rank.out = unique(rank.out), rw.evidence = rw.evidence, cs.evidence = cs.evidence, bgnet = PPI_Network))


}





### 与输入GeneList对接，生成推荐列表
omics_predict <- function(rank.predict, input_geneid, top.cutoff = 20){
  ### process # Specific diseases
  disease.gene.file <- system.file("extdata", "Disease_associated_Gene.csv", package = "dkgminer")
  disease.gene <- fread(disease.gene.file)

  diseaseid <- rank.predict$dis[1]
  spe.dis <- disease.gene[disease.gene$DiseaseId==diseaseid, ]



  ### 结合输入基因信息
  omics.rank.df<-rank.predict[rank.predict$Genes %in% input_geneid, ]

  omics.rank.df$is.dis<-rep("No",dim(omics.rank.df)[1])
  omics.rank.df[omics.rank.df$Genes %in% unique(spe.dis$GeneId),]$is.dis<-"Yes"
  omics.rank.df<-merge(omics.rank.df, unique(spe.dis[,c("GeneId", "score")]), by.x="Genes", by.y="GeneId", all.x=TRUE)

  omics.rank.df$gene.type<-rep("unsure",dim(omics.rank.df)[1])
  if( length(which( (omics.rank.df$score>0.3 & !is.na(omics.rank.df$score)) == "TRUE" )) > 0 ){
    omics.rank.df[omics.rank.df$score>0.3 & !is.na(omics.rank.df$score),]$gene.type<-"Known.dis"
  }



  ### 已知的疾病基因
  omics.rank.df1<-omics.rank.df[omics.rank.df$gene.type=="Known.dis", ]
  omics.rank.df1<-omics.rank.df1[order(omics.rank.df1$score,decreasing = TRUE),]

  omics.rank.df.k<-spe.dis[spe.dis$GeneId %in% omics.rank.df1$Genes,]
  omics.rank.df.k<-omics.rank.df.k[,c("DiseaseId","GeneId","GeneSymbol","score","pmid", "source")]
  # write.csv(omics.rank.df.k,paste("result",result.file,paste("1-Known.dis.in.omics.csv",sep = ""),sep = "/"),row.names = FALSE)


  ### 预测出的疾病基因
  omics.rank.df2<-omics.rank.df[omics.rank.df$gene.type!="Known.dis",]
  omics.rank.df2<-omics.rank.df2[order(omics.rank.df2$Score, decreasing = TRUE),]
  omics.rank.df2$Gene.sym <- mapIds(org.Hs.eg.db,
                                    keys=omics.rank.df2$Genes,
                                    keytype="ENTREZID",
                                    column="SYMBOL",
                                    multiVals="first")

  omics.rank.df.p <- na.omit(unique(omics.rank.df2[, c("Genes","Gene.sym","Score","sup.algo")]))
  if(nrow(omics.rank.df.p) > top.cutoff){ omics.rank.df.p <- omics.rank.df.p[1:top.cutoff, ] }

  omics.rank.df.p$Score <- signif(omics.rank.df.p$Score, 3)
  # write.csv(omics.rank.df.p,paste("result",result.file,paste("2-omics.rank.df.csv",sep = ""),sep = "/"),row.names = FALSE)


  return( list(known.disgene=omics.rank.df.k, predict.disgene=omics.rank.df.p) )
}
