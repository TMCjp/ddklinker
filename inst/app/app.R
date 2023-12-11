### R package app -- 0910
library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinycssloaders)
library(shinyalert)
library(DT)
library(plyr)
library(dplyr)
library(visNetwork)
library(data.table)
library(org.Hs.eg.db)
library(igraph)
library(BioNet)

#library(dkgminer)

bgnet <- 'HPRD'

### selected disease
dismat <- fread(system.file("extdata", "do_onto.csv", package = "ddklinker"))
vector_of_search_possibilities <- c('', dismat$Term)

gnid <- c('Gene symbol', 'Entrez gene ID')


### process # Specific diseases
# disease.gene.file <- system.file("extdata", "Disease_associated_Gene.csv", package = "ddklinker")
# disease.gene <- fread(disease.gene.file)
load(system.file("extdata", "DAG.Rdata", package = "ddklinker"))


### gene information
geneinfo <- fread(system.file("extdata", "geneinfo.txt", package = "ddklinker"))


### ontology information
load(system.file("extdata", "ontoinfo.rdata", package = "ddklinker"))

### graph information
load(system.file("extdata", "HPRD_network.Rdata", package = "ddklinker"))
g <- PPI_Network

### term2gene info
load(system.file("extdata", "term2gene.rdata", package = "ddklinker"))

### disease info map
load(system.file("extdata", "do2umls.rdata", package = "ddklinker"))
load(system.file("extdata", "umls2ge.rdata", package = "ddklinker"))


### network method
# netmeth <- c("RWR diffusion algorithm","Deepwalk embedding algorithm","Node2vec embedding algorithm","LINE embedding algorithm","SDNE embedding algorithm","HOPE embedding algorithm","Laplacian embedding algorithm")
netmeth <- c("RWR diffusion algorithm","Deepwalk embedding algorithm","Node2vec embedding algorithm","LINE embedding algorithm","HOPE embedding algorithm","Laplacian embedding algorithm")

### mean roc of method
mean_roc <- fread(system.file("extdata", "all_mean_roc_complete_website.csv", package = "ddklinker"))





### 0912 add disgenet link
# 添加 target="_blank" 在新标签页打开 link
createLink <- function(links, umls = 'C0002395', id, nm) {
  if(links == 'DisGeNet'){
    lk <- paste0('<a href="https://www.disgenet.org/browser/0/1/1/', umls, '/0/25/geneid__', id, '-source__ALL/" target="_blank">', nm, '</a>')
  }

  if(links == 'NCBI'){
    lk <- paste0('<a href="https://www.ncbi.nlm.nih.gov/gene/', id, '" target="_blank">', nm, '</a>')
  }

  if(links == 'pubmed'){
    lk <- paste0('<a href="https://pubmed.ncbi.nlm.nih.gov/', id, '" target="_blank">', nm, '</a>')
  }

  return(lk)
}



### 0304 add custom layout
customlay <- function(graph, conorder) {
  ig <- igraph::graph_from_data_frame(graph$x$edges[, c("from",
                                                        "to")], directed = TRUE, vertices = graph$x$nodes[,
                                                                                                          c("id", setdiff(names(graph$x$nodes), "id"))])

  coord <- igraph::layout_as_star(graph = ig, order = conorder)


  graph$x$nodes$x <- coord[, 1]
  graph$x$nodes$y <- coord[, 2]
  to <- c(-1, 1)
  from <- range(graph$x$nodes$x, na.rm = TRUE, finite = TRUE)
  if (length(unique(from)) > 1) {
    graph$x$nodes$x <- (graph$x$nodes$x - from[1])/diff(from) *
      diff(to) + to[1]
  }
  from <- range(graph$x$nodes$y, na.rm = TRUE, finite = TRUE)
  if (length(unique(from)) > 1) {
    graph$x$nodes$y <- (graph$x$nodes$y - from[1])/diff(from) *
      diff(to) + to[1]
  }
  # https://search.r-project.org/CRAN/refmans/visNetwork/html/visIgraphLayout.html
  graph$x$igraphlayout <- list(type = "square")
  # graph %>% visNodes(physics = FALSE) %>% visEdges(smooth = FALSE) %>%
  #   visPhysics(stabilization = FALSE)
  graph
}



### 1204 点击网络节点展示共有条目
comfunc <- function(x1, x2, diseaseid, path.gene, path.enrich, GO.gene, GO.enrich, pheno.gene, pheno.enrich, drug.gene, drug.enrich){
  # x1代表预测基因，x2代表已知疾病基因
  # path.gene/GO.gene 等通过加载rdata引入
  pair.term<-data.frame()

  y.path<-intersect(path.gene[(x1==path.gene$gene),]$pathway,path.gene[(x2==path.gene$gene),]$pathway)
  n.path<-length(y.path)
  if(n.path>0){
    temp<-data.frame(Genes=rep(x1,n.path),known.dis.genes=rep(x2,n.path))
    temp$term.type<-rep("pathway",n.path)
    temp$term<-y.path
    temp<-merge(temp,path.enrich[,c("ID","GeneRatio","p.adjust")],by.x="term",by.y="ID",all.x=TRUE)
    temp$is.dis<-rep("N",n.path)
    if(dim(temp[!is.na(temp$GeneRatio),])[1]>0){
      temp[!is.na(temp$GeneRatio),]$is.dis<-"Y"
    }
    temp<-temp[order(-temp$GeneRatio,temp$p.adjust),]
    temp$Rank<-c(1:dim(temp)[1])
    temp$term.e<-rep("",dim(temp)[1])
    # temp$method<-rep(method,dim(temp)[1])
    temp$dis<-rep(diseaseid,dim(temp)[1])
    # temp$n.term<-rep(n.path,dim(temp)[1])

    temp<-temp[,c("dis","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]
    # temp<-temp[,c("dis","method","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]

    pair.term<-rbind(temp,pair.term)
  }
  #GO
  # print(paste("GO term",(paste(j,N.Term,sep="/")),sep=":"))

  y.GO<-intersect(GO.gene[(x1==GO.gene$gene),]$pathway,GO.gene[(x2==GO.gene$gene),]$pathway)
  n.GO<-length(y.GO)
  if(n.GO>0){
    temp<-data.frame(Genes=rep(x1,n.GO),known.dis.genes=rep(x2,n.GO))
    temp$term.type<-rep("GO",n.GO)
    temp$term<-y.GO
    temp<-merge(temp,GO.enrich[,c("ID","GeneRatio","p.adjust")],by.x="term",by.y="ID",all.x=TRUE)
    temp$is.dis<-rep("N",n.GO)
    if(dim(temp[!is.na(temp$GeneRatio),])[1]>0){
      temp[!is.na(temp$GeneRatio),]$is.dis<-"Y"
    }
    temp<-temp[order(-temp$GeneRatio,temp$p.adjust),]
    temp$Rank<-c(1:dim(temp)[1])
    temp$term.e<-rep("",dim(temp)[1])
    # temp$method<-rep(method,dim(temp)[1])
    temp$dis<-rep(diseaseid,dim(temp)[1])
    # temp$n.term<-rep(n.GO,dim(temp)[1])

    temp<-temp[,c("dis","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]
    # temp<-temp[,c("dis","method","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]

    pair.term<-rbind(temp,pair.term)
  }

  #PHENO
  # print(paste("phenotype term",(paste(j,N.Term,sep="/")),sep=":"))
  y.pheno<-intersect(pheno.gene[(x1==pheno.gene$gene),]$pathway,pheno.gene[(x2==pheno.gene$gene),]$pathway)
  n.pheno<-length(y.pheno)
  if(n.pheno>0){
    temp<-data.frame(Genes=rep(x1,n.pheno),known.dis.genes=rep(x2,n.pheno))
    temp$term.type<-rep("pheno",n.pheno)
    temp$term<-y.pheno
    temp$is.dis<-rep(0,n.pheno)
    if(dim(temp[temp$term %in% pheno.enrich$ID,])[1]>0){
      temp[temp$term %in% pheno.enrich$ID,]$is.dis<-1
    }
    temp$GeneRatio<-rep("",n.pheno)
    temp$p.adjust<-rep("",n.pheno)

    temp<-temp[order(-(temp$is.dis),temp$term),]

    if(dim(temp[temp$is.dis==1,])[1]>0){
      temp[temp$is.dis==1,]$is.dis<-"Y"
    }
    if(dim(temp[temp$is.dis==0,])[1]>0){
      temp[temp$is.dis==0,]$is.dis<-"N"
    }

    temp$Rank<-c(1:dim(temp)[1])
    temp$term.e<-rep("",dim(temp)[1])
    # temp$method<-rep(method,dim(temp)[1])
    temp$dis<-rep(diseaseid,dim(temp)[1])
    # temp$n.term<-rep(n.pheno,dim(temp)[1])

    temp<-temp[,c("dis","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]
    # temp<-temp[,c("dis","method","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]

    pair.term<-rbind(temp,pair.term)
  }


  #drug
  # print (paste("drug term",(paste(j,N.Term,sep="/")),sep=":"))
  y.drug<-intersect(drug.gene[(x1==drug.gene$gene),]$pathway,drug.gene[(x2==drug.gene$gene),]$pathway)
  n.drug<-length(y.drug)
  if(n.drug>0){
    temp<-data.frame(Genes=rep(x1,n.drug),known.dis.genes=rep(x2,n.drug))
    temp$term.type<-rep("drug",n.drug)
    temp$term<-y.drug
    temp$is.dis<-rep(0,n.drug)
    if(dim(temp[temp$term %in% drug.enrich$ID,])[1]>0){
      temp[temp$term %in% drug.enrich$ID,]$is.dis<-1
    }
    temp$GeneRatio<-rep("",n.drug)
    temp$p.adjust<-rep("",n.drug)

    temp<-temp[order(-temp$is.dis,temp$term),]
    if(dim(temp[temp$is.dis==1,])[1]>0){
      temp[temp$is.dis==1,]$is.dis<-"Y"
    }
    if(dim(temp[temp$is.dis==0,])[1]>0){
      temp[temp$is.dis==0,]$is.dis<-"N"
    }

    temp$Rank<-c(1:dim(temp)[1])
    temp$term.e<-rep("",dim(temp)[1])
    # temp$method<-rep(method,dim(temp)[1])
    temp$dis<-rep(diseaseid,dim(temp)[1])
    # temp$n.term<-rep(n.drug,dim(temp)[1])

    temp<-temp[,c("dis","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]
    # temp<-temp[,c("dis","method","Genes","known.dis.genes","term.type","term","term.e","is.dis" ,"Rank" ,"GeneRatio","p.adjust")]

    pair.term<-rbind(temp,pair.term)

  }


  return(pair.term)
}







ui <- tagList(

  useShinyjs(),

  tags$head(
    tags$title("DDK-Linker"),
    tags$link(rel = "icon",
              type = "image/png",
              href = "fvicon.png")
  ),

  tags$style("html,body{background-color: white;overflow-x:auto;}.container{width: 100%;margin: 0 auto;padding: 0;} #myimg{width:30%;} @media screen and (min-width: 20px){.container{width: 1260px;}}"),

  tags$div(HTML(
    '<div class="container" style="width:1120px;">


<div class="row" id="top">
	<div class="col-md-3">
		<a href="http://ddklinker.ncpsb.org.cn/" target="_blank">
        <!--<img src="ddklogo.png" style="width:120%;height:80px;margin-top:5px">-->
        <img src="ddklogo.png" style="width:320px;height:80px;margin-top:5px">
        </a>
	</div>
	<!-- <div class="col-md-1"></div> -->
    <div class="col-md-9">
		<h3 class="descSite"> A user-friendly website to identify and annotate <span style="color:#984807">Disease signals</span></h3>
		<h3 class="descSite"> by linking high-throughput omics Datasets and disease Knowledge <h3>
    </div>

</div>

<style>
    .descSite{
        text-align: center;
        font-weight: bolder;
        line-height: 10px;
        color: #428bca;
    }
</style>




		</div>'
  )),

  # function page
  tags$div(class="container",

           fluidPage(
             theme = "bootstrap.css",
             includeCSS("www/styles.css"),

             navbarPage(
               # "DKGminer - Disease knowledge graph miner",
               '',
               id = "main_navbar",



               ### page 1
               tabPanel(
                 "Home",

                 dashboardPage(
                   dashboardHeader(disable = TRUE),
                   dashboardSidebar(disable = TRUE),

                   dashboardBody(

                     fluidRow(
                       box(title = 'Explore your dataset', status = "primary", solidHeader = TRUE, width = 6, height = 760,
                           collapsible = FALSE, collapsed = FALSE,


			   selectizeInput("doterm", label = tags$div(HTML("<h5><b>1. Input disease name (<a href='https://disease-ontology.org/' target='_blank'>Disease ontology</a> disease name recommended)</b> <a href='http://ddklinker.ncpsb.org.cn/help.html#Disease%20name' target='_blank'><span class='glyphicon glyphicon-question-sign'></span></a></h5>")),
                                          choices = vector_of_search_possibilities,   # c('hepatocellular carcinoma', ''),
                                          selected = NULL
                           ),

                           selectInput("geneid",
                                       label = tags$div(HTML("<h5><b>2. Select the type of gene identifier</b> <a href='http://ddklinker.ncpsb.org.cn/help.html#Select gene' target='_blank'><span class='glyphicon glyphicon-question-sign'></span></a></h5>")),
                                       choices = gnid, selectize = TRUE, selected = NULL),


                           fluidRow(
                            column(6, HTML("<h5><b>3. Input gene list for this disease</b> <a href='http://ddklinker.ncpsb.org.cn/help.html#Input gene list' target='_blank'><span class='glyphicon glyphicon-question-sign'></span></a></h5>"))
                           ),

                           textAreaInput("genels",
                                         #label = tags$div(HTML("<h5><b>3. Input gene list for this disease</b> <a href='http://ddklinker.ncpsb.org.cn/help.html#Input gene list' target='_blank'><span class='glyphicon glyphicon-question-sign'></span></a></h5>")),
                                         label = '',
                                         value = NULL, height = "150px",
                                         placeholder = "Only support human genes\nOne gene each line, or genes separated by comma"),


                           radioButtons("conf",
                                        tags$div(HTML("<h5><b>4. Set DisGeNET score for seed genes</b> <a href='http://ddklinker.ncpsb.org.cn/help.html#Set a cut off' target='_blank'><span class='glyphicon glyphicon-question-sign'></span></a></h5>")),
                                        inline = T, selected = "med",
                                        c(
                                          # "Low(no filter)" = "low",
                                          "Medium(>0.3, recommended)" = "med",
                                          "High(>0.6)" = "hig")),

                           selectizeInput("meth",
                                              label = tags$div(HTML("<h5><b>5. Select one or more linkage algorithms that you can run simultaneously</b></h5><h5 style='white-space: pre;'>                                                         <a href='http://ddklinker.ncpsb.org.cn/help.html#Network-based%20approaches' target='_blank'><u>(Performance / case study of these algorithms <sup><span class='glyphicon glyphicon-question-sign'></span></sup>)</u></a> </h5>")),
                                              choices = netmeth, multiple = TRUE),

                           # br(),
                           fluidRow(
                             column(9, ''),
                             HTML("<h5 style='display: inline-block;font-weight: bold'>(Examples:</h5>"),
                             actionButton("load1", HTML("<u>#1</u>"), style='background-color: #ffffff;color: #3c8dbc;border-color: #ffffff;padding: 10px 5px;font-size: 14px;font-weight: 900;'),
                             actionButton("load2", HTML("<u>#2</u>"), style='background-color: #ffffff;color: #3c8dbc;border-color: #ffffff;padding: 10px 0px;font-size: 14px;font-weight: 900;'),
                             HTML("<h5 style='display: inline-block;font-weight: bold'>)</h5>")
			   ),

                           fluidRow(
                             column(8, ''),
                             column(2, actionButton("reset", "Reset", style='background-color: #3C8DBC;width: 82px;'), class = "text-center"),
                            #  column(2, actionButton("load1", "Example 1", style='background-color: #3C8DBC;'), class = "text-center"),
                            #  column(2, actionButton("load2", "Example 2", style='background-color: #3C8DBC;'), class = "text-center"),
                             column(2, actionButton("explore", "Explore", style='background-color: #3C8DBC;'), class = "text-center")
                             # column(3, '')
			  )
                       ),

                       box(title = 'Brief Introduction', status = "primary", solidHeader = TRUE, width = 6, height = 760,
                           collapsible = FALSE, collapsed = FALSE,


                           div(style = "font-size: 14px;color: black;",
                               HTML(
                                 '<p><a href="http://ddklinker.ncpsb.org.cn/" target="_blank"><b><span style="color:black"><u>DDK-Linker</u></span></b></a> is a user-friendly website for Linking high-throughput Disease omics Datasets and disease Knowledge based on <a href="http://ddklinker.ncpsb.org.cn/help.html#Network-based%20approaches" target="_blank"><b><span style="color:black"><u>seven well-recognized algorithms</u></span></b></a><sup><a href="http://ddklinker.ncpsb.org.cn/help.html#Network-based%20approaches" target="_blank"><span class="glyphicon glyphicon-question-sign" style="font-size: 11px;"></span></a></sup>, to facilitate the exploration of<span style="color:#984807"> <b>disease signals</b></span><sup><a href="http://ddklinker.ncpsb.org.cn/help.html#disease%20signal" target="_blank"><span class="glyphicon glyphicon-question-sign" style="font-size: 11px;"></span></a></sup> in high-throughput screenings.</p>
                                <br>
				<img src="disease.signals.png" height="300" width="540" align="center">
				<br>
				<br>
				<p><span style="color:#984807"> <b>Disease signals</b></span><sup><a href="http://ddklinker.ncpsb.org.cn/help.html#disease%20signal" target="_blank"><span class="glyphicon glyphicon-question-sign" style="font-size: 11px;"></span></a></sup> are the most promising, biologically plausible and disease-associated genes in high-throughput disease -omics experimental data (<span style="font-size: 13px;"><a href="https://pubmed.ncbi.nlm.nih.gov/35197084/"><i>Genome Biol.</i> 2022</a>, <a href="https://pubmed.ncbi.nlm.nih.gov/31564248/"><i>Genome Med.</i> 2019</a></span>). High-throughput experimental approaches such as gene expression microarrays or quantitative proteomics allow researchers to measure <b>thousands of </b>genes, generating <b>hundreds of</b> <span style="color:#984807"><b>genes of interest</b></span> (<span style="color:#984807"> <b>GOIs</b></span><sup><a href="http://ddklinker.ncpsb.org.cn/help.html#GOIs" target="_blank"><span class="glyphicon glyphicon-question-sign" style="font-size: 11px;"></span></a></sup>, such as differentially expressed genes/proteins). However, researchers usually select only <b>several</b> genes (<span style="color:#984807"> <b>disease signals</b></span>) for further disease mechanism and drug discovery. Identifying these most valuable and informative disease signals from the massive -omics data is still a great challenge for the current biomedical community.

				Here, <a href="http://ddklinker.ncpsb.org.cn/" target="_blank"><b><span style="color:black"><u>DDK-Linker</u></span></b></a> attempts to alleviate this challenge by linking the -omics data to known disease knowledge. DDK-linker is designed to be a routine bioinformatics tool for the disease -omics experiments.<a href="http://ddklinker.ncpsb.org.cn/help.html" target="_blank"> <b>More...</b></a></p>'
                               )
                           )
                       )
                     ),



                     ## 第三行排布一个
                     fluidRow(
                       tabBox(
                         width = 12,
                         tabPanel(
                           "More Info",
                           p(h5("DDK-Linker", style='display: inline-block;font-weight: bold'),
                             h5("was developed at Beijing Proteome Research Center (BPRC).", style='display: inline-block;')
                           ),
                           p(h5("Browser compatibility:", style='display: inline-block;font-weight: bold'),
                             h5("We support the latest versions of the" , style='display: inline-block;'),
                             a(href="https://chorm.com.cn/", 'Chrome', target="_blank", style="font-size: 13.28px;"),
                             h5(',', style='display: inline-block;'),
                             a(href="https://www.mozilla.org/en-US/firefox/new/", 'Firefox', target="_blank", style="font-size: 13.28px;"),
                             h5(',', style='display: inline-block;'),
                             a(href="https://www.opera.com/", 'Opera', target="_blank", style="font-size: 13.28px;"),
                             h5(',', style='display: inline-block;'),
                             a(href="https://www.apple.com/safari/", 'Safari', target="_blank", style="font-size: 13.28px;"),
                             h5(',', style='display: inline-block;'),
                             a(href="https://www.microsoft.com/en-us/download/internet-explorer.aspx", 'Internet Explorer', target="_blank", style="font-size: 13.28px;"),
                             h5("(version 10 and 11)." , style='display: inline-block;')
                           ),
                           p(h5("Citation:", style="display: inline-block;font-weight: bold;"),
                             h5("Diao L, Jiang P, Kong X,", style="display: inline-block;font-size: 13.28px;"),
                             h5("et al. ,", style='display: inline-block;font-style:italic;font-size: 13.28px;'),
                             h5("DDK-Linker: a user-friendly website to identify and annotate disease signals from high-throughput -omics datasets in the context of ", style="display: inline-block;font-size: 13.28px;"),
                             h5("disease knowledge graph.", style="display: inline-block;font-size: 13.28px;"),
                             h5("(submitted to", style="display: inline-block;font-size: 13.28px;"),
                             h5("Nucleic Acids Research", style='display: inline-block;font-style:italic;font-weight: bold;font-size: 13.28px;'),
                             h5(")", style='display: inline-block;font-size: 13.28px;'))
                         ),
                         tabPanel(
                           "Acknowledgements",
                           tags$div(
                             HTML('<p><b>We gratefully acknowledge <a href="https://www.disgenet.org/">DisGeNET</a>, <a href="http://www.hprd.org/">HPRD</a> for their generously sharing of their datasets and algorithm.</b><br><br>
                                        <b>This work is funded by the:</b><br>
										National Natural Science Foundation of China [32088101, 32271518,31871341], National key Research and Development Program of China [2016YFC0901905, 2020YFE0202200, 2021YFA1301603] and Program of International S&T Cooperation (2014DFB30020).<br>
									<img src="ISTC.jpg" style="width: 80px;">&nbsp;&nbsp;&nbsp;<img src="NNSFC.jpg" style="width: 200px;">&nbsp;&nbsp;&nbsp;<img src="973.jpg" style="width: 60px;"></p>')
                           )),
                         tabPanel(
                           "Website History",
                           tags$div(
                             HTML(
                               '<div style="margin:10px 10px 5px 10px;">
                           <p><span style="color:red">#2022.12.20</span> Proposal of DDK-Liner 1.0 was submitted to <a href="https://academic.oup.com/nar"><i>Nucleic Acids Research WSI</i></a>.</p>
					  				       <p><span style="color:red">#2022.12.01</span> 130 successful analyses have been conducted by outside users.</p>
					  				       <p><span style="color:red">#2022.06.01</span> DDK-Linker was officially online for outside users.</p>
					  				       <p><span style="color:red">#2022.04.01</span> The underlying disease knowledge graph was expanded to cover 2000 human diseases.</p>
					  				       <p><span style="color:red">#2022.03.01</span> Network embedding analysis function was added</p>
					  				       <p><span style="color:red">#2022.01.01</span> DDK-Linker starts internal testing, used to find the candidate disease gene by network fusion analysis.</p>
					  			         </div>'
                             )
                           )
                         )
                       )

                     )

                   )

                 )
               ),



               ### page 2
               tabPanel(
                 "My analysis",
                 # icon("home",lib = "glyphicon"),
                 fluidRow(actionButton("back", HTML("<u><< Back</u>"), style='background-color: #ffffff;color: #3c8dbc;border-color: #ffffff;padding: 5px 10px;font-size: 18px;font-weight: 900;')),
                 withSpinner(uiOutput("respage"))
               ),

               ### page 5
               tabPanel(
                 "Case study",

                 #  ### add to homepage
                 #  HTML(
                 #    '<style>
                 #              .asite {text-decoration: underline;}
                 #            </style>

                 #            <a class="asite" href="http://ddklinker.ncpsb.org.cn/" target="_blank" style="font-size: 16px;color: #4384BD;font-weight: bold;"> << Back to Homepage </a>
                 #           '
                 #  ),
                 htmlOutput("case_frame")
               ),


               ### page 4
               tabPanel(
                 "Documentation",
                 htmlOutput("frame")
               ),


               ## page 3
               tabPanel(
                 "About",
                 htmlOutput("about_frame")
               )


               #  tabPanel(
               #    "Contact",
               #    HTML(
               #      '<style>
               #                .asite {text-decoration: underline;}
               #              </style>

               #              <a class="asite" href="http://ddklinker.ncpsb.org.cn/" target="_blank" style="font-size: 16px;color: #4384BD;font-weight: bold;"> << Back to Homepage </a>
               #             '
               #    ),
               #    htmlOutput("contact_frame")
               #  )

             )

           )
  )
)









# Define server logic required to draw a histogram
server <- function(input, output, session) {
  ### home to analysis page
  ### 初始状态各参数值
  #comment_value <- reactiveValues(disnm = '', ingns = '')
  comment_value <- reactiveValues(disnm = '', ingns = '')


  ### Event observers (三类事件：点击加载例子、用户手动输入、点击重置)
  ### 主要参考：https://stackoverflow.com/questions/45436603/shiny-initial-textareainput-value-and-reactive-after-each-button-press
  # example button
  observeEvent(input$load1,
               {
                 comment_value$disnm <- "Alzheimer's disease"
                 comment_value$ingns <- "LYN,CD69,EIF4G1,PLXNA4,SNAP29,BCR,PPP1R9B,ICA1,TXLNA,BANK1,ARHGEF12,AXIN1,INPPL1,CLIP2,CASP3,TDRKH,IKBKG,MESD,STK4,ITGB1BP2,CALCOCO1,SRPK2,DAPP1,DAB2,ZBTB16,GRAP2,SRC,SNAP23,MAP4K5,ERBIN,YES1,BACH1,SH2B3,FKBP1B,WASF1,AIFM1,MAP2K6,TRIM5,PRTFDC1,CDKN1A,PMVK,FOXO1,USO1,HEXIM1,GOPC,AIMP1,TBCB,CA13,TANK,TACC3,NFATC1,LAT2,SCAMP3,METAP1D,CBL,STX6,CRKL,DECR1,SMAD1,IRAK1,FKBP5,PTPN1,IRAK4,KIFBP,LRMP,VPS53,PLA2G4A,HSPB1,PPP1R2,NAA10,STX16,SPRY2,EGF,DCTN1,ABL1,MANF,PTPN6,FLI1,DRG2,GP6,CETN2,FGF2,LAT,PPIB,F11R,YTHDF3,MYO9B,NUB1,MGMT,GFER,FOXO3,PECAM1,CD2AP,PRDX5,NT5C3A,PRKCQ,VPS37A,PRDX3,MAX,ENO2,WWP2,CERT1,NF2,LACTB2,HCLS1,FXYD5,CASP2,LAP3,TOP2B,ANXA11,ARHGAP25,SERPINB6,PVALB,GAMT,PTPRJ,ARHGAP1,TBL1X,AKR1B1,FES,PLXNB3,BAG6,NFKBIE,SULT1A1,COMT,CDC27,ILKAP,STX8,RRM2B,HTRA2,AKT1S1,VASH1,TRAF2,BIRC2,EIF4B,IQGAP2,FADD,HMOX2,RP2,RPS6KB1,IMPA1,METAP2,LGALS8,WAS,CRADD,DCTN2,DFFA,SELP,SIRT2,CD63,STAMBP,TYMP,DAG1,DIABLO,STXBP3,P4HB,CD40LG,NUDT5,PRKRA,FHIT,BGN,TP53,PSME1,KYAT1,WASF3,CLEC1B,USP8,MIF,IRF9,PARK7,EDAR,DGKZ,BTC,SCARF1,MVK,ERP44,DNAJB1,LIFR,ARSB,MAGED1,TMSB10,ANXA4,QDPR,PRDX1,AHCY,PRKAB1,PAG1,GSAP,CCT5,STIP1,VAMP5,HDGF,KYNU,INPP1,GLB1,ACAA1,MCFD2,PAK4,ENAH,SH2D1A,FKBP7,PLXDC1,TXNDC5,BID,MAEA,CXCL1,F2R,CCL5,ITGB1BP1,EGLN1,CD164,TIGAR,ATP6V1D,AIF1,RASSF2,EIF5A,PEBP1,DPP7,PPM1B,LGMN,GALNT2,FKBP4,CD84,PIK3AP1,PRDX6,CNTN5,GP1BA,ITGA6,NAMPT,ATG4A,PFDN2,CALR,DDX58,CD40,SUMF2,BLMH,CAMKK1,KLK4,CXCL5,TCL1A,PFKM,FGR,TPP1,STC1,NUCB2,LAMA4,TRIM21,ING1,PTX3,PPP3R1,ABHD14B,EGFR,MMP7,MEP1B,ITGB7,LRP1,AOC3,CD8A,ATP6V1F,NADK,PTP4A1,IL1B,HSPB6,SKAP1,HPGDS,SPINK4,CNPY2,CD46,IGSF3,PLAU,DKK4,CRELD2,FAP,KLK14,CD97,RET,FETUB,TNFSF13B,PAPPA,CSF1,THOP1,ITGB1,KRT19,GLO1,SOD2,SERPINE1,MMP3,ALDH1A1,FGF5,TNFAIP8,PDP1,SMOC1,GUSB,DPP10,AGRP,PSIP1,ITGB2,FUT8,DEFB4A,MASP1,SIRT5,CX3CL1,APBB1IP,ENTPD2,DCTPP1,CSNK1D,SDC4,AARSD1,CRHBP,ITGA11,PHOSPHO1,TNC,CFC1,CNTN2,SDC1,HBEGF,TGFA,CTRC,WNT9A,CCL17,C1QA,BRK1,NCS1,ANXA1,LTA4H,CDHR5,NRTN,SEPTIN9,DPEP1,CTF1,CCL11,GALNT10,ROBO2,FAM3B,CHL1,DDC,CCL2,IL13RA1,FGFBP1,PCSK9,OSMR,IL7,ALCAM,CDON,SIGLEC7,PDGFA,IFNLR1,CDH17,ACP5,DPP4,EIF4EBP1,PARP1,IL1R2,TNFSF10,NCF2,TNFSF14,FLT1,XCL1,TNFRSF14,SCG2,CHIT1,PXN,CES2,VCAM1,BAMBI,SOD1,CCN1,NBN,VAT1,EZR,ERBB2,ACTN4,COCH,FUS,DCN,ESAM,NFATC3,APEX1,NELL1,TNFRSF10B,PRSS2,ERBB3,METAP1,PPY,CBLN4,UMOD,HNMT,MMP1,CNDP1,SNCG,CTSD,SCLY,PDGFRA,MICB,ADM,OMG,TIMP4,CANT1,ANGPTL4,AREG,NOMO1,CDH5,S100A11,FAS,TNFRSF10A,CPM,VEGFD,AOC1,FLT3,FABP9,MANSC1,PLA2G10,GFRA1,PDGFB,EPHA10,IGFBP3,IGFBP2,TGFBR3,FBP1,CLSTN2,FGF19,PAM,CLSPN,TFRC,ULBP2,TNR,LYPD1,CNTN1,PREB,ZBTB17"

                 updateTextAreaInput(session, "doterm", value = comment_value$disnm)
                 updateTextAreaInput(session, "genels", value = comment_value$ingns)

                 updateTextAreaInput(session, "geneid", value = 'Gene symbol')
                 updateTextAreaInput(session, "conf", value = 'med')

		 #updateSelectInput(session, "meth", choices = c("RWR (power:0.8615)", "Deepwalk (power:0.8276)", "Node2Vec (power:0.7134)", "LINE (power:0.5837)", "HOPE (power:0.8212)", "Laplacian (power:0.7432)"))
                 #updateTextAreaInput(session, "meth", value = c("RWR (power:0.8615)", "Deepwalk (power:0.8276)"))
               }
  )

  observeEvent(input$load2,
               {
                 comment_value$disnm <- "Alzheimer's disease"
                 comment_value$ingns <- "AGRN,CR1,NCK2,BIN1,INPP5D,CLNK,TINP1,HAVCR2,HLA-DRB1,TREM2,CD2AP,TMEM106B,ZCWPW1,NYAP1,EPHA1-AS1,CLU,SHARPIN,USP6NL,ECHDC3,CCDC6,MADD,SPI1,MS4A4A,PICALM,SORL1,FERMT2,RIN3,ADAM10,APH1B,SCIMP,RABEP1,GRN,ABI3,TSPOAP1-AS1,ACE,ABCA7,APOE,NTN5,CD33,LILRB2,CASS4,APP"

                 updateTextAreaInput(session, "doterm", value = comment_value$disnm)
                 updateTextAreaInput(session, "genels", value = comment_value$ingns)

                 updateTextAreaInput(session, "geneid", value = 'Gene symbol')
                 updateTextAreaInput(session, "conf", value = 'med')

		 #updateSelectInput(session, "meth", choices = c("RWR (power:0.8615)", "Deepwalk (power:0.8276)", "Node2Vec (power:0.7134)", "LINE (power:0.5837)", "HOPE (power:0.8212)", "Laplacian (power:0.7432)"))
                 #updateTextAreaInput(session, "meth", value = c("RWR (power:0.8615)", "Deepwalk (power:0.8276)"))
               }
  )


  # 在上面load例子的时候dignm会被赋值，相当于是用户输入了，所以下面的功能会启动，就造成了meth部分会刷新，用户的meth例子不能加载上。231209进行了修改，上面加载例子的时候不包括meth，meth在下面刷新后，再默认选择前两个即可
  # user input (选择项或输入文本来自用户的输入)
  observeEvent(input$doterm, {
    comment_value$disnm <- input$doterm
    updateTextAreaInput(session, "doterm", value = comment_value$disnm)

    nn <- match(input$doterm, dismat$Term)
    if(!is.na(nn)){
      diseaseid <- dismat$do_id[nn]

      cn <- match(diseaseid, mean_roc$dis.Name)
      if(!is.na(cn)){
        myChoices <- unlist(lapply(1:6, function(x) paste0(colnames(mean_roc)[x+1], ' (power:', sprintf("%0.4f", as.numeric(mean_roc[cn, -1])[x]), ')')))
      } else {
        myChoices <- c('RWR','Deepwalk','Node2Vec','LINE','HOPE','Laplacian')
      }

      updateSelectInput(session, "meth", choices = myChoices)
      updateTextAreaInput(session, "meth", value = myChoices[1:2])   #在meth刷新后再默认选择前两个meth
    }

  })


  observeEvent(input$genels, {
    comment_value$ingns <- input$genels
    updateTextAreaInput(session, "genels", value = comment_value$ingns)
  })

  # reset button
  # https://rdrr.io/cran/shinyjs/man/reset.html
  observeEvent(input$reset, {
    # reset("parameter")
    comment_value$disnm <- ''
    comment_value$ingns <- ''

    updateTextAreaInput(session, "doterm", value = comment_value$disnm)
    updateTextAreaInput(session, "genels", value = comment_value$ingns)
    updateTextAreaInput(session, "geneid", value = 'Gene symbol')
    updateTextAreaInput(session, "conf", value = 'med')
    updateSelectInput(session, "meth", choices = netmeth)
    updateTextAreaInput(session, "meth", value = '')
  })



  hideTab(inputId = "main_navbar", target = "My analysis")


  ### back to home
  observeEvent(input$back, {
    updateTabsetPanel(session, "main_navbar", selected = "Home")
  })



  ### return args
  eargs <- eventReactive(input$explore, {

    ### submit info
    runjs("let doterm = $('#doterm').val();let genels=$('#genels').val(); let params={doterm, genels};$.ajax({url:'http://111.198.139.89/spaceIPRecords/save',type:'POST',dataType:'json',data:{'websiteTag':'63b595bb0b02bc00373f7128','params':JSON.stringify(params)}})")



    nn <- match(input$doterm, dismat$Term)
    if(!is.na(nn)){
      diseaseid <- dismat$do_id[nn]

      ### 2006 ava vector
      if(dismat$available[nn] == 'F'){
        alerthtml <- HTML(
          paste0(
            '<p style="line-height: 22px; font-size: 15px;">Sorry, the current version DDK-linker could not continue the analysis for your designated disease <span style="color:#984807"><b>"', input$doterm, '"</b></span>, because there are no high-confidence known disease genes for this disease in <a href="https://www.disgenet.org/" target="_blank">DisGeNET</a>.</p>
              <br>
              <p style="line-height: 22px; font-size: 15px;">To ensure the accuracy of the analysis results, we only link your GOI to high confidence disease genes in DisGeNET. If you want to know more details, please contact us by email <a href="mailto:lidong.bprc@163.com">lidong.bprc@foxmail.com</a>. </p>
          <br>
          <p style="line-height: 22px; font-size: 15px;">Here are <a href="http://ddklinker.ncpsb.org.cn/Human.diseases.in.DDK-Linker.xlsx" download="Human.diseases.in.DDK-Linker.xlsx">2006 diseases</a> with high confidence disease genes. Please refer to <a href="http://ddklinker.ncpsb.org.cn/about_iframe.html#dis_info" target="_blank">Human diseases in DDK-Linker</a> for more details.</p>'
          )
        )

        shinyalert("", tagList( div(style = "text-align: left;", alerthtml) ),
                   type = "", size = "s", html = TRUE,
                   showConfirmButton = F, immediate = T)
      }
    }





    genetype <- input$geneid

    genelists <- input$genels
    if(length(grep('\n', genelists, perl = T)) > 0){
      genelists <- strsplit(genelists, split = '\n', perl = T)[[1]]
    }
    if(length(grep(',', genelists, fixed = T)) > 0){
      genelists <- strsplit(genelists, split = ',', fixed = T)[[1]]
    }

    genelists <- gsub(' ', '', genelists, fixed = T)
    genelists <- toupper(genelists)

    if(genetype == "Gene symbol"){
      input_geneid <- mapIds(org.Hs.eg.db,
                             keys=genelists,
                             keytype="SYMBOL",
                             column="ENTREZID",
                             multiVals="first")
    } else {
      input_geneid <- genelists
    }


    ###
    #meth <- input$meth
    #methmap <- cbind(c("RWR diffusion algorithm","Deepwalk embedding algorithm","Node2vec embedding algorithm","LINE embedding algorithm","SDNE embedding algorithm","HOPE embedding algorithm","Laplacian embedding algorithm"), c("RW","deepwalk","node2vec","line","sdne","HOPE","Laplacian"))
    #meth <- methmap[match(meth, methmap), 2]


    ### 230530
    meth <- input$meth
    meth <- unlist( lapply(meth, function(x) strsplit(x, split = ' ', fixed = T)[[1]][1]) )
    methmap <- cbind(c("RWR","Deepwalk","Node2Vec","LINE","HOPE","Laplacian"), c("RW","deepwalk","node2vec","line","HOPE","Laplacian"))
    meth <- methmap[match(meth, methmap), 2]


    ###
    ### predict disease genes
    # rfd <- paste0('D://Project/R_project/dkgminer/inst/app/www/', gsub(':', '.', diseaseid, fixed = T), '.rdata')
    if(input$conf == 'low'){cutd = 'cutoff.0.0'}
    if(input$conf == 'med'){cutd = 'cutoff.0.3';cutsel = 0.3}
    if(input$conf == 'hig'){cutd = 'cutoff.0.6';cutsel = 0.6}

    rfd <- url(paste0('http://ddklinker.ncpsb.org.cn/netdata/', cutd, '/', gsub(':', '.', diseaseid, fixed = T), '.', bgnet, '.rdata'))
    load(rfd)
    close(rfd)


    ### output args
    list(diseaseid=diseaseid, meth=meth, input_geneid=input_geneid,
         evi=candicate.evidence.df, rank=candicate.rank.df,
         cutsel=cutsel,
         go_gs=GO.gene, path_gs=path.gene, drug_gs=drug.gene, pheno_gs=pheno.gene,
         go=GO.enrich, path=path.enrich, drug=drug.enrich, pheno=pheno.enrich)
  })






  ### run DKG analysis
  observeEvent(input$explore, {
    if(input$doterm == '' | input$genels == ''){
      # showNotification("Please complete the parameters", type = "error", duration = 10, action = a(href = "javascript:location.reload();", "Reload page"))
      shinyalert("", "Warning! You did not submit any gene list. \nPlease submit again.", type = "warning", size = "s", showConfirmButton = F, immediate = F)
    } else {
      # runjs("var today = new Date(); alert(today);")
      shinyalert("", "Performing lingkage analysis,this analysis takes about 10 seconds", type = "success", size = "xs", showConfirmButton = F, immediate = F, timer = 30000)
      shinyalert("", "Analysis Error, please check the parameters and restart", type = "error", size = "xs", showConfirmButton = F, immediate = F)
    }



    ### import args
    diseaseid <- eargs()$diseaseid
    meth <- eargs()$meth
    input_geneid <- eargs()$input_geneid
    cutsel <- eargs()$cutsel

    candicate.evidence.df <- eargs()$evi
    candicate.rank.df <- eargs()$rank



    ##### known disease genes
    ### 23.2.6 cutoff=0.1
    spe.dis <- disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score >= 0.1)), ]
    nn <- unlist(lapply(input_geneid, function(x) which(spe.dis$GeneId == x)))
    disgen <- spe.dis[nn, ]

    ### 23.2.15
    intdg3 <- length(unique(disgen[which(disgen$score >= 0.3), ]$GeneId))
    intdg6 <- length(unique(disgen[which(disgen$score >= 0.6), ]$GeneId))


    ### gene2dis
    disgen <- merge(do2umls, disgen, by.x='DOID', by.y='DiseaseId', all.y = T)
    disgen <- merge(umls2ge, disgen, by = c('diseaseId', 'GeneSymbol', 'score'), all = F)
    disgen <- disgen[order(disgen$score, decreasing = T), ]

    # if(input$conf == 'low'){cutv = 0}
    # if(input$conf == 'med'){cutv = 0.3}
    # if(input$conf == 'hig'){cutv = 0.6}
    # 新增候选基因二次划分（disgenet中打分小于0.1的保留候选基因，其他归为已知疾病基因）
    # cutv <- 0.1
    # nn <- which(disgen$score >= cutv)
    # disgen <- disgen[nn, ]

    cutgs <- unique(disgen$GeneId)


    disgen$score <- createLink('DisGeNet', disgen$diseaseId, disgen$GeneId, disgen$score)
    disgen$GeneSymbol <- createLink('NCBI', id = disgen$GeneId, nm = disgen$GeneSymbol)
    #disgen$pmid <- createLink('pubmed', disgen$pmid, disgen$pmid)


    ##### known disease genes
    # disgen <- read.table('D:/Project/R_project/RpackTest/netdat/dkgdata/1-Known.dis.in.omics.csv', sep = ',', header = T)
    disgen <- disgen[, c('GeneSymbol', 'score', 'pmid')]

    gs <- unique(disgen$GeneSymbol)
    disgenDT <- as.data.frame(matrix(NA, nrow = length(gs), ncol = 3))

    if(length(gs) > 0){
      for(i in 1:length(gs)){
        nn <- match(gs[i], disgen$GeneSymbol)

        disgenDT[i, 1] <- disgen[nn, 1]
        disgenDT[i, 2] <- disgen[nn, 2]
        # disgenDT[i, 3] <- disgen[nn, 3]

        pmn <- disgen$pmid[grep(gs[i], disgen$GeneSymbol)]
        if(length(pmn) == 1 && length(which(is.na(pmn))) > 0){
          pmn <- ''
        } else {
          pmn <- unique(na.omit(pmn))
          pmn <- createLink('pubmed', id = pmn, nm = pmn)
        }

        if(length(pmn) > 3){
          #disgenDT[i, 3] <- paste0(gsub(',', ';', toString( pmn[1:3] )), ' ……')
          disgenDT[i, 3] <- gsub(';', ',', toString( pmn[1:3] ))
        } else if(length(pmn) <= 3 && length(pmn) >= 1){
          disgenDT[i, 3] <- gsub(';', ',', toString( pmn ))
        } else if(length(pmn) < 1){
          disgenDT[i, 3] <- ''
        }

      }
    }

    # colnames(disgenDT) <- c('Gene symbol', 'Score in DisGNet', 'Evidence')
    colnames(disgenDT) <- c('Gene symbol', 'DisGeNet score', 'Publications (PMID)')

    ### output known genes table
    # https://rstudio.github.io/DT/options.html，DOM Elements，f表示检索框，t表示表格, i表示汇总信息，p表示翻页选项……
    # https://github.com/rstudio/DT/issues/773，上面所述的例子如，dom = 'Bfrtip'
    output$disgen <- renderDT(disgenDT, escape = FALSE, selection = 'none',
                              options = list(pageLength = 10, dom = 'lftp', scrollX = T,
                                             initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#3c8dbc', 'color': '#fff'});","}"),
                                             # ### https://rstudio.github.io/DT/options.html
                                             # columnDefs = list(list(
                                             #     targets = 4,
                                             #     render = JS(
                                             #         "function(data) {",
                                             #         "return data.length > 10 ?",
                                             #         "'<span title=\"' + data + '\">' + data.substr(0, 10) + '...</span>' : data;",
                                             #         "}")
                                             # )),

                                             # https://github.com/rstudio/DT/issues/783
                                             headerCallback = DT::JS(
                                               "function(thead) {",
                                               "  $(thead).css('font-size', '1em');",
                                               "}")
                              ))


    ### statics info
    ### 提交基因大于0，数据库中疾病基因小于5个，不做检验只列出数字
    if(intdg3 > 0){
      stainfo3 <- 'y'


      dgenes <- length( unique( disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score >= 0.3)), ]$GeneId ) )

      if(dgenes <= 5){
        annotext3 <- 'n'

        output$statis3 <- renderText({HTML(gsub(', )', ')',
                                                paste0("<p style='font-size:13px;white-space: pre;'>        #Genes (DisGNet Score>0.3): <b>", intdg3, '</b></p>'),
                                                fixed = T)
        )})

      } else {
        annotext3 <- 'y'

        t2g <- unique(disease.gene[which(disease.gene$score >= 0.3), c("DiseaseId", "GeneId")])
        ### 讨论 如果现存算法有问题（背景基因是所有疾病基因），背景基因改成全部人类基因（来自org.HS包）
        x = clusterProfiler::enricher(input_geneid, TERM2GENE=t2g, TERM2NAME=NA, qvalueCutoff=1, pvalueCutoff=1, minGSSize=0, maxGSSize = 50000)
        tt <- x@result
        tt <- tt[which(tt$ID == diseaseid), ]

        grl <- as.numeric(strsplit(tt$GeneRatio, '/', fixed = T)[[1]])
        gr <- grl[1]/grl[2]
        mapnum1 <- grl[1]
        brl <- as.numeric(strsplit(tt$BgRatio, '/', fixed = T)[[1]])
        br <- brl[1]/brl[2]
        ER <- round(gr/br, digits = 3)
        pv1 <- signif(tt$pvalue, 3)

        if(ER>1 && pv1<0.05){sign1 = 'significantly enriched'; ervalue1 = ER}
        if(ER>1 && pv1>=0.05){sign1 = 'enriched'; ervalue1 = ER}
        if(ER<=1){sign1 = ''; ervalue1 = '-'}

        if(pv1>=0.05){
          sign1_html <- paste0('<span style="color:#984807">', sign1, '</span>')
        }
        if(pv1<0.05){
          sign1_html <- paste0('<span style="color:red">', sign1, '</span>')
        }

        output$statis3 <- renderText({HTML(
          paste0("<p style='font-size:13px;white-space: pre;'>        #Genes (DisGNet Score>0.3): <b>", mapnum1, "</b>, ", sign1_html, " (Enrichment ratio=", ervalue1, ', p-value=', pv1, '*)</p>')
        )})
      }
    } else {
      stainfo3 <- 'n'
      annotext3 <- 'n'
    }



    if(intdg6 > 0){
      stainfo6 <- 'y'


      dgenes <- length( unique( disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score >= 0.6)), ]$GeneId ) )

      if(dgenes <= 5){
        annotext6 <- 'n'

        output$statis6 <- renderText({HTML(gsub(', )', ')',
                                                paste0("<p style='font-size:13px;white-space: pre;'>        #Genes (DisGNet Score>0.6): <b>", intdg6, '</b></p>'),
                                                fixed = T)
        )})
      } else {
        annotext6 <- 'y'

        t2g <- unique(disease.gene[which(disease.gene$score >= 0.6), c("DiseaseId", "GeneId")])
        ### 讨论 如果现存算法有问题（背景基因是所有疾病基因），背景基因改成全部人类基因（来自org.HS包）
        x = clusterProfiler::enricher(input_geneid, TERM2GENE=t2g, TERM2NAME=NA, qvalueCutoff=1, pvalueCutoff=1, minGSSize=0, maxGSSize = 50000)
        tt <- x@result
        tt <- tt[which(tt$ID == diseaseid), ]

        grl <- as.numeric(strsplit(tt$GeneRatio, '/', fixed = T)[[1]])
        gr <- grl[1]/grl[2]
        mapnum2 <- grl[1]
        brl <- as.numeric(strsplit(tt$BgRatio, '/', fixed = T)[[1]])
        br <- brl[1]/brl[2]
        ER <- round(gr/br, digits = 3)
        pv2 <- signif(tt$pvalue, 3)

        if(ER>1 && pv2<0.05){sign2 = 'significantly enriched'; ervalue2 = ER}
        if(ER>1 && pv2>=0.05){sign2 = 'enriched'; ervalue2 = ER}
        if(ER<=1){sign2 = ''; ervalue2 = '-'}

        if(pv2>=0.05){
          sign2_html <- paste0('<span style="color:#984807">', sign2, '</span>')
        }
        if(pv2<0.05){
          sign2_html <- paste0('<span style="color:red">', sign2, '</span>')
        }

        output$statis6 <- renderText({HTML(paste0("<p style='font-size:13px;white-space: pre;'>        #Genes (DisGNet Score>0.6): <b>", mapnum2, "</b>, ", sign2_html, " (Enrichment ratio=", ervalue2, ', p-value=', pv2, '*)</p>'))})
      }
    } else {
      stainfo6 <- 'n'
      annotext6 <- 'n'
    }









    ##### predicte genes
    nn <- unlist(lapply(meth, function(x) which(candicate.rank.df$method == x)))
    candicate.rank.df <- candicate.rank.df[nn, ]

    nn <- na.omit(unlist(lapply(input_geneid, function(x) which(candicate.rank.df$Genes == x))))
    dagres_all <- candicate.rank.df[nn, ]

    dagres <- matrix(NA, nrow = length(unique(dagres_all$Genes)), ncol = 3)
    dagres[, 1] <- unique(dagres_all$Genes)
    for(r in 1:length(unique(dagres_all$Genes))){
      cadg <- dagres[r, 1]
      dagres[r, 2] <- min(as.numeric(dagres_all$Rank[which(dagres_all$Genes == cadg)]))

      ### 0306, mean top2000 rank
      rknum <- as.numeric(dagres_all$Rank[which(dagres_all$Genes == cadg)])
      rn <- which(rknum <= 2000)
      if(length(rn) > 0){rknum <- rknum[rn]}
      ###

      dagres[r, 3] <- mean(rknum)
    }
    colnames(dagres) <- c('Genes', 'Top.rank', 'Mean.rank')

    dagres <- as.data.frame(dagres)
    dagres$Top.rank <- as.numeric(dagres$Top.rank)

    ### 0306, top rank 1000
    tn <- which(dagres$Top.rank <= 1000)
    if(length(tn) > 0){dagres <- dagres[tn, ]}
    ###

    dagres$Mean.rank <- sprintf("%0.1f", as.numeric(dagres$Mean.rank))
    dagres <- dagres[order(dagres$Top.rank, decreasing = F), ]

    dagres_all <- dagres

    dagres <- cbind(dagres, tx=rep('Detail', nrow(dagres)))
    dagres <- unique(dagres)


    if(nrow(dagres) > 0){
      dagres$Symbol <- mapIds(org.Hs.eg.db,
                              keys=dagres$Genes,
                              keytype="ENTREZID",
                              column="SYMBOL",
                              multiVals="first")
      dagres <- dagres[, c("Genes","Symbol","Top.rank","Mean.rank","tx")]

      cn <- na.omit(match(cutgs, dagres$Genes))
      #if(length(cn) > 0){dagres <- dagres[-cn, ]}
      if(length(cn) > 0){
        dagres <- dagres[-cn, ]
        dagres_all <- dagres_all[-cn, ]
      }
    }

    if(nrow(dagres) > 0){
      ### output result number
      # https://itecnote.com/tecnote/r-outputting-multiple-lines-of-text-with-rendertext-in-r-shiny/
      # output$num <- renderText({ paste0('Candidate disease associated genes in submitted genes (', nrow(dagres), ')') })
      # output$num <- renderText({ paste0('2 Known disease associated genes (', nrow(disgenDT), ') in submitted genes') })
      output$num <- renderText({
        paste0('<h3 style="font-weight: bold;display: inline-block;">2 Known disease-associated genes ', "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Confidence Score' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 12px;'></span></a></sup>", ' (', nrow(disgenDT), ') in submitted genes</h3>')
      })


      ### 23.2.14   Asterisks(*) denote that the marked candidate genes are low confidence(<0.1) disease-associated genes in DisGeNET
      low.ges <- unique(disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score < 0.1)), ]$GeneId)
      nn <- na.omit(match(low.ges, dagres$Genes))
      if(length(nn) > 0){dagres$Symbol[nn] <- paste0('*', dagres$Symbol[nn])}


      rownames(dagres) <- 1:nrow(dagres)
      colnames(dagres) <- c("id", "Gene", "Top.rank", "Mean.rank", "Details")

      ### 1220 temp
      # dagres$Gene <- createLink('NCBI', dagres$id, dagres$Gene)
      dagres <- dagres[, -1]


      ### 0309
      dagres$Top.rank <- paste0(signif(as.numeric(dagres$Top.rank)/20000, digits = 2)*100, '%')
      dagres$Mean.rank <- paste0(signif(as.numeric(dagres$Mean.rank)/20000, digits = 2)*100, '%')


      output$dagtable <- renderDataTable(server = FALSE,{ datatable(dagres, escape = FALSE, selection = list(mode = "single", selected = 1),
                                                                    # # https://blog.csdn.net/weixin_42780204/article/details/118636045
                                                                    extensions = 'Buttons',
                                                                    options = list(pageLength = 10, dom = 'Bfrtip',   #buttons = c('csv', 'copy'),
                                                                                   buttons =list(list(extend = "csv",                                                                                            filename = "Table_of_candidate_disease-associated_genes",
                                                                                                      text = '<span class="glyphicon glyphicon-download-alt"></span> CSV',
                                                                                                      exportOptions = list(columns = c(1,2,3),modifier = list(page = "all"))
                                                                                   )),
                                                                                   columnDefs = list(list(className = 'dt-center', targets = c(2,3))),
                                                                                   initComplete = JS("function(settings, json) {","$(this.api().table().header()).css({'background-color': '#3c8dbc', 'color': '#fff'});","}")
                                                                    )) %>%
          formatStyle(4, cursor = 'pointer') })


      ### compute network files
      visnetfile <- reactive({
        # input是独立的参数，表示shiny页面选择的输入
        stb = input$dagtable_rows_selected

        if(!is.null(stb)){
          sg <- dagres_all$Genes[stb]
        } else {
          stb <- 1
          sg <- dagres_all$Genes[stb]
        }


        net.df <- candicate.evidence.df[which(candicate.evidence.df$Genes == sg), 2:4]
        gn <- unlist(lapply(net.df$method, function(x) length(grep(x, meth, ignore.case = T)) > 0))
        net.df <- net.df[gn, ]

	### del >top5
        nmath <- unlist(lapply(unique(net.df$method), function(x) if(length(which(net.df$method == x)) > 5){which(net.df$method == x)[1:5]} else {which(net.df$method == x)} ))
        net.df <- net.df[nmath, ]


        ngid <- unique(c(net.df$Genes, net.df$known.dis.genes))
        ngsy <- mapIds(org.Hs.eg.db,
                       keys=ngid,
                       keytype="ENTREZID",
                       column="SYMBOL",
                       multiVals="first")
        net.node <- cbind(id=1:length(ngid), label=ngsy)
        net.node <- as.data.frame(net.node)


        #diao0106 add xgmml download
        net.node.output<-as.data.frame(ngsy)
        net.node.output$geneid<-row.names(net.node.output)
        names(net.node.output)<-c("Genesym","Geneid")
        net.node.output<-net.node.output[,c(2,1)]
        net.node.output$Gene.type <- c('candi.gene', rep('known.dis', nrow(net.node) - 1))

        net.df.output<-as.data.frame(net.df[,c(2,3,1)])
        names(net.df.output)[1]<-"candidate.genes"
        net.df.output$method.group<-rep("other",dim(net.df.output)[1])
        if(dim(net.df.output[net.df.output$method=="rw",])[1]>0){
          net.df.output[net.df.output$method=="rw",]$method.group<-net.df.output[net.df.output$method=="rw",]$method
        }

        net <- igraph::graph_from_data_frame(net.df.output, vertices=net.node.output)
        # saveNetwork(network = net, name="network", file =paste(path.temp,meth,".",sg,sep = "") , type="XGMML")
        # end


        net.edge <- net.df[, c(2,3,1)]
        colnames(net.edge) <- c('from', 'to', 'title')
        net.edge$from <- mapvalues(net.edge$from, from = ngid, to = 1:length(ngid), warn_missing = F)
        net.edge$to <- mapvalues(net.edge$to, from = ngid, to = 1:length(ngid), warn_missing = F)


        ### 0304 ring node rank
        starid <- unique(net.edge$to)
        starconnet <- unlist(lapply(starid, function(x) length(net.edge$title[which(net.edge$to == x)])))
        starmat <- cbind(id=starid, connet=starconnet)

        conorder <- merge(net.node[, 1:2], starmat, by='id', all.x = T)
        conorder[1,3] <- 10
        conorder$connet <- as.numeric(conorder$connet)

        conorder <- conorder[order(conorder$label), ]
        # conorder <- conorder[order(conorder$connet, decreasing = TRUE), ]
        conorder <- conorder$id





        if(length(which(net.edge$title == 'rw')) > 0){
          net.edge_eb <- net.edge[-which(net.edge$title == 'rw'), ]
        } else {
          net.edge_eb <- net.edge
        }
        if(nrow(net.edge_eb) > 0){
          net.edge_cb <- unique(net.edge_eb[, 1:2])
          methtx <- unlist(lapply(net.edge_cb$to, function(x) toString(net.edge_eb$title[which(net.edge_eb$to == x)]) ))
          net.edge_cb <- cbind(net.edge_cb, title=methtx)

          net.edge <- rbind(net.edge[which(net.edge$title == 'rw'), ], net.edge_cb)
        }


        ### nodes info
        # net.node$pro.type <- NA
        net.node$pro.type <- c('candi.gene', rep('known.dis', nrow(net.node) - 1))

        net.node$color.background[which(net.node$pro.type == 'candi.gene')] <- '#984807'   #'orange'
        net.node$color.background[grep('known.dis', net.node$pro.type)] <- '#006600'   #'#FFFFFF'

        # https://github.com/datastorm-open/visNetwork/issues/51
        net.node$color.border[which(net.node$pro.type == 'candi.gene')] <- '#984807'
        net.node$color.border[grep('known.dis', net.node$pro.type)] <- '#006600'

        net.node$shape[grep('known.dis', net.node$pro.type)] <- 'dot'
        net.node$shape[grep('candi.gene', net.node$pro.type)] <- 'circle'

        # https://rdrr.io/cran/shp2graph/man/Degree.list.html
        # nodes$value <- shp2graph::Degree.list(nodes[, 1:2], edges[, 1:2], Directed = F)[[1]]

        ### 通过比例缩放设定每个节点的大小： https://github.com/datastorm-open/visNetwork/issues/49
        #net.node$value[grep('known.dis', net.node$pro.type)] <- 1.2
        #net.node$value[grep('candi.gene', net.node$pro.type)] <- 1.5
	net.node$value[grep('candi.gene', net.node$pro.type)] <- 2
    	for(x in 1:nrow(starmat)){
		net.node$value[which(net.node$id == starmat[x, 1])] <- 1+0.2*as.numeric(starmat[x, 2])
    	}

        ###
        net.node$font.color[grep('known.dis', net.node$pro.type)] <- '#000000'
        net.node$font.color[grep('candi.gene', net.node$pro.type)] <- '#FFFFFF'

        net.node$font.size[grep('known.dis', net.node$pro.type)] <- 18
        net.node$font.size[grep('candi.gene', net.node$pro.type)] <- 24

        # nodes$font.size <- 8
        net.node$shadow <- FALSE



        ### edges info
        net.edge$dashes <- TRUE
        net.edge$dashes[grep('rw', net.edge$title, ignore.case = T)] <- FALSE

        net.edge$color <- '#984801'
        net.edge$color[grep('rw', net.edge$title, ignore.case = T)] <- '#18723D'

        # net.edge$width <- 3
        nmeth <- unlist(lapply(net.edge$title, function(x) stringr::str_count(x, ',')))
        net.edge$width <- 3 + 2*nmeth

        net.edge$smooth.type <- ''
        net.edge$smooth.type[grep('rw', net.edge$title, ignore.case = T)] <- "curvedCW"



        output$mt <- renderUI(
          p(h4('Graph for the linkages between ', style="display: inline-block;font-weight: bold;"),
            span(net.node$label[1], style="color: #984807;font-weight: bold;font-size: 22px;"),
            h4(' and known disease genes', style="display: inline-block;font-weight: bold;"))
        )
        output$st <- renderUI(
          p(h6("    Clicking each node to view biological annotation of the corresponding linkage.", style="white-space: pre;"),
            h6("    Known disease genes are arranged clockwise in alphabetical order of their gene symbols.", style="white-space: pre;"))
            # h6("    Know disease genes clockwise in descending order according to the number of linkage algorithms (secondary: alphabetical).", style="white-space: pre;"))
            # h5("    The thickness of the line indicates the number of algorithms that support one linkage.", style="white-space: pre;"))
        )



        ### 11.18 temp
        net.node <- na.omit(net.node)

        ### network output
        output$dhtm <- downloadHandler(
          filename = function(){
            #paste0("Network-", Sys.Date(), ".html")
            paste0("Linkages_between_",net.node$label[1],"_and_known_disease_genes",".html")

          },
          content = function(file){
            visNetwork(net.node, net.edge) %>% visSave(file, background = "white")
          }
        )

        output$dxgm <- downloadHandler(
          filename = function(){
            paste0("Linkages_between_",net.node$label[1],"_and_known_disease_genes",".xgmml")
          },
          content = function(file){
            saveNetwork(network = net, name="network", file, type="XGMML")
          }
        )

        # https://9to5answer.com/r-shiny-download-existing-file
        output$dxml <- downloadHandler(
          filename = function(){
            'net.styles.xml'
          },
          content = function(file){
            file.copy("http://ddklinker.ncpsb.org.cn/DDK-Linker.net.styles.xml", file)
          }
        )

        vg <- visNetwork(net.node, net.edge, width = "100%") %>%   #, main = mt, submain = st) %>%
          visNodes( scaling = list(label = list(enabled = T, min = 16, max = 18)) )

        vg <- vg %>% customlay(conorder)

        vg %>%
          # visNodes( scaling = list(label = list(enabled = T, min = 8, max = 15)) ) %>%
          # visEdges(smooth = list(type = "curvedCW")) %>%
          # visExport(type = "png", label = 'PNG', style= "color: #ffffff;background: #3C90BC;border:none;") %>%
          visOptions(manipulation = FALSE, highlightNearest = FALSE, nodesIdSelection = list(enabled = FALSE)) %>%
          # visLegend(addNodes = lnodes, addEdges = ledges, position = "right", useGroups = FALSE, ncol = 1) %>%
          visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = TRUE) %>%
          visPhysics(maxVelocity = 0.1) %>%
          visEvents(selectNode = "function(nodes) {Shiny.onInputChange('current_node_id', nodes);}")
      })


      output$visnet <- renderVisNetwork({
        visnetfile()
      })


      ### analysis success info
      # shinyalert("", "The analysis is complete. See the results", type = "success", size = "s", showConfirmButton = T, confirmButtonText = "OK", immediate = T)
      shinyalert("", "The analysis is complete. See the results", type = "success", size = "s", timer = 1, immediate = T)








      ### Generate UI (直接通过server生成前端页面 result页面，防止插件嵌套冲突)
      output$respage <- renderUI({
        fluidPage(

          fluidRow(
            column(12,

                   HTML(
                     paste0('<h2 style="display: inline-block;font-weight: bold;">Explore disease-associated signals ',
                            "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#output' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 12px;'></span></a></sup>",
                            ' in submitted genes for ',
                            '<span style="color: #984807;font-weight: bold;font-size: 32px;">', input$doterm, '</span></h2>')
                   ),
                   hr(style="height:2px;border-width:0;background-color:gray"),

                   h3('1 Candidate disease-associated genes in submitted genes', style="font-weight: bold;"),

                   #  HTML(
                   #    paste0(
                   #      '<h5 style="white-space: pre;display: inline-block;">    All the known disease-associated genes',
                   #      "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Network diffusion strategy' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 8px;'></span></a></sup>",
                   #      ' are omitted from the prediction results.', '</h5>'
                   #    )
                   #  ),
                   HTML(
                     paste0(
                       '<h5 style="white-space: pre;display: inline-block;">    All your selected seed genes (DisGeNET score >', cutsel, ') are omitted from the table of candidate disease-associated gene',
                       "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Network diffusion strategy' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 8px;'></span></a></sup>",
                       ' .', '</h5>'
                     )
                   ),

                   h5('    Click on "Details" in the left table to view the corresponding linkages in the figure of the right panel.', style="white-space: pre;"),
                   h5('    Mean.rank is the mean rank across the ranks of selected algorithms (normalize). Top.rank is the maximum scaled rank assigned to each gene across selected algorithms (normalize).', style="white-space: pre;")
            )
          ),

          useShinyalert(force = TRUE),

          fluidRow(
            column(width = 4,
                   div(
                     ### TODO:balance
                     p(' ', style="line-height:5"),
                     h4('Table of candidate disease-associated genes', style="text-align: left;font-weight: bold;white-space: pre;color: #984807;"),
                     p(),
                     DTOutput("dagtable"),
                     HTML("<p style='font-size:13px'>Asterisks(*) denote that the marked candidate genes are low confidence(<0.1) disease-associated genes in DisGeNET, and they are not used as seed genes for linkage analysis(<a href='https://www.frontiersin.org/articles/10.3389/fgene.2021.779186/full' target='_blank'>as Du et al. defined</a>)</p>")
                   )
            ),
            column(width = 7, offset = 1,
                   # https://community.rstudio.com/t/fluidrow-and-column-add-border-to-the-respective-block/13187/2
                   style = "border: 1px solid black;",
                   # br(),
                   uiOutput("mt"),
                   uiOutput("st"),
                   fluidRow(
                     column(6),
                     #    column(10),
                     #column(1, downloadButton("dsvg", "SVG", style='color: black; background-color: #FFFFFF; padding:6px; font-size:80%', icon = icon("download")), class = "text-center"),
                     column(1, downloadButton("dhtm", "HTML", style='color: black; background-color: #FFFFFF; padding:6px; font-size:80%', icon = icon("download")), class = "text-center"),
                     column(2, downloadButton("dxgm", "XGMML", style='color: black; background-color: #FFFFFF; padding:6px; font-size:80%', icon = icon("download")), class = "text-center"),
                     column(1, downloadButton("dxml", "Cytoscape VizMap", style='color: black; background-color: #FFFFFF; padding:6px; font-size:80%', icon = icon("download")), class = "text-center")
                   ),

                   visNetworkOutput("visnet",height = "420px") %>% withSpinner(type = 5),
                   img(src = "legend.png", height = 80, width = 600, align = "center")
            )
          ),

          hr(style="height:2px;border-width:0;background-color:gray"),

          fluidRow(
            tags$div(htmlOutput("num")),
            if(stainfo6 == 'y' | stainfo3 == 'y'){HTML("<p style='line-height: 12px;font-size:13px;white-space: pre;'> <b>Statistics</b>:<br></p>")},
            if(stainfo6 == 'y'){htmlOutput("statis6")},
            if(stainfo3 == 'y'){htmlOutput("statis3")},
            if(annotext3 == 'y' | annotext6 == 'y'){HTML("<p style='line-height: 12px;font-size:13px;white-space: pre;'>        * One-tailed hypergeometric distribution test</p>")},
            br(),
            DTOutput("disgen"),
            br()
          )
        )
      })
    } else {
      ### output result number
      output$num <- renderText({
        paste0('<h3 style="font-weight: bold;display: inline-block;">Known disease-associated genes ', "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Confidence Score' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 12px;'></span></a></sup>", ' (', nrow(disgenDT), ') in submitted genes</h3>')
      })


      alerthtml <- HTML(
        '<p style="line-height: 22px; font-size: 14px;">All genes of interest you submitted are known disease-associated genes in <a href="https://www.disgenet.org/" target="_blank">DisGeNET</a>. Here we will only present the details of these genes without performing linkage analysis.</p>'
      )
      shinyalert("", tagList( div(style = "text-align: left;", alerthtml) ), html = TRUE, type = "success", size = "s", timer = 5000, immediate = T)


      output$respage <- renderUI({

        fluidPage(
          fluidRow(
            column(12,

                   HTML(
                     paste0('<h2 style="display: inline-block;font-weight: bold;">Explore disease-associated signals ',
                            "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#output' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 12px;'></span></a></sup>",
                            ' in submitted genes for ',
                            '<span style="color: #984807;font-weight: bold;font-size: 32px;">', input$doterm, '</span></h2>')
                   ),
                   hr(style="height:2px;border-width:0;background-color:gray")
            )
          ),


          useShinyalert(force = TRUE),

          # titlePanel( textOutput("num") ),
          fluidRow(
            # h3(textOutput("num"), style="font-weight: bold;"),
            tags$div(htmlOutput("num")),
            # p("5 known disease associated genes (DisGNet confidence score >0.3) in submitted genes"),

            # https://stackoverflow.com/questions/24049159/change-the-color-and-font-of-text-in-shiny-app
            # tags$head(tags$style("#num{color: black;font-size: 16px;font-weight: bold;}")),
            if(stainfo6 == 'y' | stainfo3 == 'y'){HTML("<p style='line-height: 12px;font-size:13px;white-space: pre;'> <b>Statistics</b>:<br></p>")},
            if(stainfo6 == 'y'){htmlOutput("statis6")},
            if(stainfo3 == 'y'){htmlOutput("statis3")},
            if(annotext3 == 'y' | annotext6 == 'y'){HTML("<p style='line-height: 12px;font-size:13px;white-space: pre;'>        * One-tailed hypergeometric distribution test</p>")},

            br(),
            DTOutput("disgen"),
            br()
          )
        )
      })
    }






    ### to show result tabpanel
    showTab(inputId = "main_navbar", target = "My analysis")
    # updateTabsetPanel(session, "main_navbar", selected = "Home")
    updateTabsetPanel(session, "main_navbar", selected = "My analysis")
    # updateTabsetPanel(session, "main_navbar", selected = "Analysis")
  })




  ### click network nodes
  # https://deanattali.com/blog/shinyalert-package/
  observeEvent(input$current_node_id, {
    shinyalert("", "Please wait 3 seconds", type = "success", size = "s", showConfirmButton = F, immediate = F)

    diseaseid <- eargs()$diseaseid
    meth <- eargs()$meth
    input_geneid <- eargs()$input_geneid

    candicate.evidence.df <- eargs()$evi
    candicate.rank.df <- eargs()$rank

    GO.enrich <- eargs()$go
    path.enrich <- eargs()$path
    drug.enrich <- eargs()$drug
    pheno.enrich <- eargs()$pheno
    GO.gene <- eargs()$go_gs
    path.gene <- eargs()$path_gs
    drug.gene <- eargs()$drug_gs
    pheno.gene <- eargs()$pheno_gs


    spe.dis <- disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score >= 0.1)), ]
    nn <- unlist(lapply(input_geneid, function(x) which(spe.dis$GeneId == x)))
    disgen <- spe.dis[nn, ]

    cutv <- 0.1
    nn <- which(disgen$score >= cutv)
    disgen <- disgen[nn, ]

    cutgs <- unique(disgen$GeneId)


    ##### predicte genes
    nn <- unlist(lapply(meth, function(x) which(candicate.rank.df$method == x)))
    candicate.rank.df <- candicate.rank.df[nn, ]

    nn <- na.omit(unlist(lapply(input_geneid, function(x) which(candicate.rank.df$Genes == x))))
    dagres_all <- candicate.rank.df[nn, ]

    dagres <- matrix(NA, nrow = length(unique(dagres_all$Genes)), ncol = 3)
    dagres[, 1] <- unique(dagres_all$Genes)
    for(r in 1:length(unique(dagres_all$Genes))){
      cadg <- dagres[r, 1]
      dagres[r, 2] <- min(as.numeric(dagres_all$Rank[which(dagres_all$Genes == cadg)]))

      ### 0306, mean top2000 rank
      rknum <- as.numeric(dagres_all$Rank[which(dagres_all$Genes == cadg)])
      rn <- which(rknum <= 2000)
      if(length(rn) > 0){rknum <- rknum[rn]}
      ###

      dagres[r, 3] <- mean(rknum)
    }
    colnames(dagres) <- c('Genes', 'Top.rank', 'Mean.rank')

    dagres <- as.data.frame(dagres)
    dagres$Top.rank <- as.numeric(dagres$Top.rank)

    ### 0306, top rank 1000
    tn <- which(dagres$Top.rank <= 1000)
    if(length(tn) > 0){dagres <- dagres[tn, ]}
    ###

    dagres$Mean.rank <- as.numeric(dagres$Mean.rank)
    dagres_all <- dagres[order(dagres$Top.rank, decreasing = F), ]


    cn <- na.omit(match(cutgs, dagres_all$Genes))
    if(length(cn) > 0){
      dagres_all <- dagres_all[-cn, ]
    }





    stb = input$dagtable_rows_selected

    if(!is.null(stb)){
      sg <- dagres_all$Genes[stb]
    } else {
      stb <- 1
      sg <- dagres_all$Genes[stb]
    }

    net.df <- candicate.evidence.df[which(candicate.evidence.df$Genes == sg), 2:4]
    gn <- unlist(lapply(net.df$method, function(x) length(grep(x, meth, ignore.case = T)) > 0))
    net.df <- net.df[gn, ]

    ngid <- unique(c(net.df$Genes, net.df$known.dis.genes))

    net.node <- cbind(id=1:length(ngid), label=ngid)
    net.node <- as.data.frame(net.node)


    s <- input$current_node_id$nodes[[1]]

    if(s > 1){
      dg <- net.node$label[which(net.node$id == s)]

      # showNotification("Please wait 3 seconds while searching", type = "message", duration = 3)
      functab <- comfunc(x1 = sg, x2 = dg, diseaseid, path.gene, path.enrich, GO.gene, GO.enrich, pheno.gene, pheno.enrich, drug.gene, drug.enrich)

      ### shortest path
      short.path<-shortest_paths(g, from = V(g)[name == sg], to =  V(g)[name=dg], output = 'both')
      short.path.v <-short.path$vpath
      if( length(short.path.v[[1]]) > 0 ){
        short.path.v1 <- sapply(short.path.v, as_ids)
        short.path.v1 <- short.path.v1[, 1]

        shortG <- unlist( lapply(short.path.v1, function(x) geneinfo$`Approved symbol`[which(geneinfo$`NCBI Gene ID` == x)]) )
      }


      ### disgenet
      nn <- match(input$doterm, dismat$Term)
      if(!is.na(nn)){ diseaseid <- dismat$do_id[nn] }
      spe.dis <- disease.gene[intersect(which(disease.gene$DiseaseId==diseaseid), which(disease.gene$score >= 0.3)), ]

      ### gene2dis
      disgen <- merge(do2umls, spe.dis, by.x='DOID', by.y='DiseaseId', all.y = T)
      disgen <- merge(umls2ge, disgen, by = c('diseaseId', 'GeneSymbol', 'score'), all = F)

      dagsinfo <- unique(disgen[, c('diseaseId', 'GeneId', 'GeneSymbol', 'score')])
      dagsinfo$link <- createLink('DisGeNet', dagsinfo$diseaseId, dagsinfo$GeneId, dagsinfo$score)

      dagsinfo <- dagsinfo$link[which(dagsinfo$GeneId == dg)]


      # sgsym <- mapIds(org.Hs.eg.db, keys=sg, keytype="ENTREZID", column="SYMBOL", multiVals="first")
      # dgsym <- mapIds(org.Hs.eg.db, keys=dg, keytype="ENTREZID", column="SYMBOL", multiVals="first")
      sgsym <- geneinfo$`Approved symbol`[which(geneinfo$`NCBI Gene ID` == sg)]
      dgsym <- geneinfo$`Approved symbol`[which(geneinfo$`NCBI Gene ID` == dg)]

      sgname <- geneinfo$`Approved name`[which(geneinfo$`NCBI Gene ID` == sg)]
      dgname <- geneinfo$`Approved name`[which(geneinfo$`NCBI Gene ID` == dg)]


      ### ontology link
      nn <- grep('^path', functab$term, perl = T)
      npath <- gsub('path:', '', functab$term[nn], fixed = T)
      if(length(npath) > 0){
        keggname <- ontoinfo$name[match(npath, ontoinfo$term)]
        keggls <- lapply(1:length(npath), function(x)
          div(
            h5('        ', style="display: inline-block; white-space: pre;"),
            # a(href=paste0("https://www.kegg.jp/pathway/", npath[x], '+', sg, '+', dg), img(src = "pathway.png", width = "25px", height = "20px"), target="_blank"),
            a(href=paste0("https://www.kegg.jp/pathway/", npath[x], '+', sg, '+', dg), paste0(npath[x], ': ', keggname[x]), target="_blank", style="font-size: 12px;color: #1e90ff;")
            # a(href=paste0("https://www.kegg.jp/pathway/", npath[x], '+', sg, '+', dg), img(src = paste0("kegg.path.png/", npath[x], ".png"), width = "32px", height = "28px", style="border: 2px solid #3C8DBC;"), target="_blank")

          )
        )
      }


      nn <- grep('^GO', functab$term, perl = T)
      ngo <- functab$term[nn]
      if(length(ngo) > 0){
        goname <- ontoinfo$name[match(ngo, ontoinfo$term)]
        gols <- lapply(1:length(ngo), function(x) a(href=paste0("http://amigo.geneontology.org/amigo/search/bioentity?q=*:*&fq=isa_partof_closure:%22", ngo[x], "%22&sfq=document_category:%22bioentity%22"), paste0(ngo[x], ': ', goname[x]), target="_blank", style="font-size: 12px;color: #1e90ff;"))
      }

      nn <- grep('^HP', functab$term, perl = T)
      nhp <- functab$term[nn]
      if(length(nhp) > 0){
        hpname <- ontoinfo$name[match(nhp, ontoinfo$term)]
        hpls <- lapply(1:length(nhp), function(x) a(href=paste0("https://hpo.jax.org/app/browse/term/", nhp[x]), paste0(nhp[x], ': ', hpname[x]), target="_blank", style="font-size: 12px;color: #1e90ff;"))
      }

      nn <- grep('^DB', functab$term, perl = T)
      ndrug <- functab$term[nn]
      if(length(ndrug) > 0){
        drugname <- ontoinfo$name[match(ndrug, ontoinfo$term)]
        drugls <- lapply(1:length(ndrug), function(x) a(href=paste0("https://go.drugbank.com/drugs/", ndrug[x]), paste0(ndrug[x], ': ', drugname[x]), target="_blank", style="font-size: 12px;color: #1e90ff;"))
      }
      ###


      ### method & sorce
      meth <- gsub('^RW$', 'rw', meth, perl = T)
      nn <- intersect(unlist(lapply(meth, function(x) which(candicate.evidence.df$method == x))), intersect(which(candicate.evidence.df$Genes == sg), which(candicate.evidence.df$known.dis.genes == dg)))
      gescore <- signif(candicate.evidence.df$Score.evi[nn], 3)
      gemeth <- gsub('^rw$', 'RWR diffusion algorithm', candicate.evidence.df$method[nn], perl = T)

      gemeth <- gsub('^deepwalk$', 'Deepwalk embedding algorithm', gemeth, perl = T)
      gemeth <- gsub('^node2vec$', 'Node2vec embedding algorithm', gemeth, perl = T)
      gemeth <- gsub('^line$', 'LINE embedding algorithm', gemeth, perl = T)
      gemeth <- gsub('^sdne$', 'SDNE embedding algorithm', gemeth, perl = T)
      gemeth <- gsub('^HOPE$', 'HOPE embedding algorithm', gemeth, perl = T)
      gemeth <- gsub('^Laplacian$', 'Laplacian embedding algorithm', gemeth, perl = T)

      strategy <- lapply(1:length(gemeth), function(x) p(style="line-height:2", h5(paste0('    ', gemeth[x], ' (', gescore[x], ');'), style="display: inline-block;white-space: pre;")))




      alerthtml <- list(
        h3(paste0('Information of linkage for ', sgsym, ' and ', dgsym),  style="background-color: #1b6fb5;color: #fff;padding: 18px 18px"),
        br(),

        h4('1 Brief information', style="color: #984807;font-weight: bold"),
        p(h5('    Candidate disease-associated gene: ', style="display: inline-block; white-space: pre;"),
          a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sgsym), paste0(sgsym, ', ', sgname), target="_blank", style="font-size: 12px;color: #1e90ff;")
        ),
        p(h5('    Known disease-associated gene: ', style="display: inline-block; white-space: pre;"),
          a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", dgsym), paste0(dgsym, ', ', dgname), target="_blank", style="font-size: 12px;color: #1e90ff;"),
          HTML(paste0('<h6 style="display: inline-block; white-space: pre;font-weight: bold"> (DisGeNET score: ', dagsinfo, ') </h6>'))
        ),
        br(),

        h4('2 Linkage algorithm', style="color: #984807;font-weight: bold"),
        strategy,
        br(),

        h4('3 Bioinformatics annotations of linkage', style="color: #984807;font-weight: bold"),
        # h5(paste0('  ', toString(functab$term)))

        h5('    (1) Shared pathways', style="white-space: pre;font-weight: bold"),
        h6('        Click the link to view the linkage for gene pairs in KEGG pathway graph', style="white-space: pre"),
        if(exists("keggls")){lapply(keggls, function(x) x)},
        if(!exists("keggls")){h5('        --', style="white-space: pre;")},

        h5('    (2) Shared gene ontology terms', style="white-space: pre;font-weight: bold"),
        if(exists("gols")){
          lapply(gols, function(x) list(h5('        ', style="display: inline-block; white-space: pre;"), x, br()))
        } else {
          h5('        --', style="white-space: pre;")
        },

        h5('    (3) Shared related drugs', style="white-space: pre;font-weight: bold"),
        if(exists("drugls")){
          lapply(drugls, function(x) list(h5('        ', style="display: inline-block; white-space: pre;"), x, br()))
        } else {
          h5('        --', style="white-space: pre;")
        },

        h5('    (4) Shared phenotype terms', style="white-space: pre;font-weight: bold"),
        if(exists("hpls")){
          if(length(hpls) > 5){
            p(
              lapply(hpls[1:5], function(x) list(h5('        ', style="display: inline-block; white-space: pre;"), x, br())),
              h5('        ……', style="white-space: pre;")
            )

          } else {
            lapply(hpls, function(x) list(h5('        ', style="display: inline-block; white-space: pre;"), x, br()))
          }
        } else {
          h5('        --', style="white-space: pre;")
        },


        if(exists("shortG")){
          p(
            h5('    (5) Protein-protein interaction network', style="white-space: pre;font-weight: bold"),
            h5('        ', style="display: inline-block; white-space: pre;"),
            h5(shortG[1], style="border: 3px solid #A52A2A;display: inline-block;"),
            if(length(shortG) > 2){
              h5(paste0('--', gsub(', ', '--', toString(shortG[2:(length(shortG)-1)]), fixed = T), '--'), style="display: inline-block;")
            } else {
              h5('--', style="display: inline-block;")
            },
            h5(shortG[length(shortG)], style="border: 3px solid #228B22;display: inline-block;")
          )
        }
      )

      # Show a modal when the button is pressed
      shinyalert("", closeOnClickOutside = TRUE, #callbackJS = 'window.scrollTo(0,0);',
                 text = tagList(
                   div(style = "text-align: left;font-family: ''Source Sans Pro',sans-serif';color: #000000;",   #font-family: 'Roboto';
                       alerthtml) ),
                 html = TRUE, type = "", size = "m", immediate = T,
                 showConfirmButton = T, confirmButtonText = 'Close', confirmButtonCol = '#1B6FB5')

    }



    ### 点击中间节点的情况
    if(s == 1){
      # shinyalert("", "Please wait 3 seconds while searching", type = "success", size = "l", showConfirmButton = F, immediate = T)
      # methnames <- gsub('^RW$', 'Random walk with restart', meth, perl = T)
      methnames <- gsub('^RW$', 'RWR diffusion algorithm', meth, perl = T)
      methnames <- gsub('^deepwalk$', 'Deepwalk embedding algorithm', methnames, perl = T)
      methnames <- gsub('^node2vec$', 'Node2vec embedding algorithm', methnames, perl = T)
      methnames <- gsub('^line$', 'LINE embedding algorithm', methnames, perl = T)
      methnames <- gsub('^sdne$', 'SDNE embedding algorithm', methnames, perl = T)
      methnames <- gsub('^HOPE$', 'HOPE embedding algorithm', methnames, perl = T)
      methnames <- gsub('^Laplacian$', 'Laplacian embedding algorithm', methnames, perl = T)

      sgsym <- geneinfo$`Approved symbol`[which(geneinfo$`NCBI Gene ID` == sg)]
      sgname <- geneinfo$`Approved name`[which(geneinfo$`NCBI Gene ID` == sg)]

      ### brief
      sn <- which(dagres_all$Genes == sg)
      nn <- nrow(dagres_all)
      sortinfo <- paste0(sn, '/', nn)
      # scoreinfo <- signif(candicate.rank.df$Score[which(candicate.rank.df$Genes == sg)], 3)

      Mdag <- dagres_all[order(dagres_all$Mean.rank, decreasing = F), ]
      sn <- which(Mdag$Genes == sg)
      nn <- nrow(Mdag)
      Msortinfo <- paste0(sn, '/', nn)





      ### method & sorce
      nn <- which(candicate.rank.df$Genes == sg)

      # rwscore <- candicate.rank.df$Score[nn]
      rwscore <- unlist( lapply(nn, function(x) signif((1-as.numeric(candicate.rank.df$Rank[x])/max(as.numeric(candicate.rank.df$Rank))), 3)) )

      rwsort <- unlist(lapply(nn, function(x) if(as.numeric(candicate.rank.df$Rank[x]) <= 2000){paste0( candicate.rank.df$Rank[x], '/', max(as.numeric(candicate.rank.df$Rank)) )} else {'-'}))



      meth <- gsub('^RW$', 'rw', meth, perl = T)
      mm <- lapply(meth, function(x) intersect(which(candicate.evidence.df$method == x), intersect(which(candicate.evidence.df$Genes == sg), which(candicate.evidence.df$Rank.evi < 6))))
      ### del >top5
      mm <- lapply(mm, function(x) if(length(x) > 5){x <- x[1:5]}else{x <- x} )

      strategy <- alist()
      for(ii in 1:length(mm)){
        nn <- mm[[ii]]
        scorename <- 'Similarity score'
        if(methnames[ii] == 'RWR diffusion algorithm'){scorename <- 'Connectivity score'}

        topgene <- candicate.evidence.df$known.dis.genes[nn]
        topgsym <- geneinfo$`Approved symbol`[match(topgene, geneinfo$`NCBI Gene ID`)]
        topscore <- signif(candicate.evidence.df$Score.evi[nn], 3)
        rwtop <- alist()
        for(x in 1:length(topgsym)){
          rwtop[[1+2*(x-1)]] <- a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", topgsym[x]), topgsym[x], target="_blank", style="font-size: 12px;color: #1e90ff;")
          if(x < length(topgsym)){rwtop[[2+2*(x-1)]] <- h5(', ', style="color: #1e90ff;display: inline-block;")}
        }

        strategy[[ii]] <- list(
          div(
            h5(paste0('    (', ii, ') Algorithm : '), style="display: inline-block;white-space: pre;"),
            h6(methnames[ii], style="display: inline-block;")
          ),
          HTML(
            paste0(
              '<h5 style="white-space: pre;display: inline-block;">        ', scorename,
              "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Network diffusion strategy' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 8px;'></span></a></sup>",
              ' : <h6 style="display: inline-block;">', rwscore[ii], ' (', rwsort[ii], ')', '</h6></h5>'
            )
          ),
          p(h5(paste0('        Top ', length(topgsym), ' known disease genes with highest connectivity to ', sgsym, ' : '), style="white-space: pre;"),
            h5('        ', style="display: inline-block; white-space: pre;"),
            rwtop
          )
        )

      }



      ###
      alerthtml <- list(
        # tags$style("html,body{background-color: white;overflow-x:auto;}.container{width: 100%;margin: 0 auto;padding: 0;} #myimg{width:30%;} @media screen and (min-width: 20px){.container{width: 1440px;}}"),

        h3(paste0('Information for candidate disease gene ', sgsym), style="background-color: #1b6fb5;color: #fff;padding: 18px 18px"),
        br(),

        h4('1 Brief information', style="color: #984807; font-weight: bold"),
        p(h5('    Candidate disease-associated gene: ', style="display: inline-block; white-space: pre;"),
          a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sgsym), paste0(sgsym, ', ', sgname), target="_blank", style="font-size: 12px;color: #1e90ff;")
        ),

        div(
          h5('    Top rank of candidate gene: ', style="display: inline-block;white-space: pre;"),
          h6(sortinfo, style="display: inline-block;"),
          h5('; Mean rank of candidate gene: ', style="display: inline-block;white-space: pre;"),
          h6(Msortinfo, style="display: inline-block;")
        ),

        # HTML(
        #   paste0(
        #     '<h5 style="white-space: pre;display: inline-block;">    Confidence score',
        #     "<sup><a href='http://ddklinker.ncpsb.org.cn/help.html#Network diffusion strategy' target='_blank'><span class='glyphicon glyphicon-question-sign' style='font-size: 8px;'></span></a></sup>",
        #     ' : <h6 style="display: inline-block;">', scoreinfo, '</h6></h5>'
        #   )
        # ),
        # h5(paste0('    Confidence score: ', scoreinfo), style="white-space: pre;"),
        div(
          h5('    Linkage algorithm: ', style="display: inline-block;white-space: pre;"),
          h6(gsub(', ', '; ', toString(gsub(' embedding algorithm', '', methnames, fixed = T)), fixed = TRUE), style="display: inline-block;")
        ),
        br(),

        h4('2 Details information of linkage algorithm', style="color: #984807; font-weight: bold"),
        strategy,
        br(),

        h4('3 Additional links', style="color: #984807; font-weight: bold"),
        p(
          h5('    NCBI Gene: ', style="white-space: pre;display: inline-block;"),
          a(href=paste0("https://www.ncbi.nlm.nih.gov/gene/", sg), sgsym, target="_blank", style="font-size: 12px;color: #1e90ff;")
        ),
        p(
          h5('    GeneCards: ', style="white-space: pre;display: inline-block;"),
          a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", sgsym), sgsym, target="_blank", style="font-size: 12px;color: #1e90ff;")
        ),
        p(
          h5('    UniProt: ', style="white-space: pre;display: inline-block;"),
          a(href=paste0("https://www.uniprot.org/uniprotkb?query=", sgsym), sgsym, target="_blank", style="font-size: 12px;color: #1e90ff;")
        )
      )

      # Show a modal when the button is pressed
      shinyalert("", closeOnClickOutside = TRUE,
                 text = tagList(
                   div(style = "text-align: left;font-family: ''Source Sans Pro',sans-serif';color: #000000;",
                       alerthtml) ),
                 html = TRUE, type = "", size = "m", immediate = T,
                 showConfirmButton = T, confirmButtonText = 'Close', confirmButtonCol = '#1B6FB5')

    }
  })







  ### https://gist.github.com/cecilialee/46e5ec16237cf17bf9a5e971edb332f7
  ### iframe help web
  output$frame <- renderUI({
    helpweb <- tags$iframe(src="help_iframe.html", style="height:12500px; width:105%", frameborder = "0", scrolling="no")
    helpweb
  })

  output$about_frame <- renderUI({
    helpweb <- tags$iframe(src="about_iframe.html", style="height:3200px; width:104%", frameborder = "0", scrolling="no")
    helpweb
  })

  output$contact_frame <- renderUI({
    helpweb <- tags$iframe(src="contact.html", style="height:520px; width:100%", frameborder = "0", scrolling="no")
    helpweb
  })

  output$case_frame <- renderUI({
    helpweb <- tags$iframe(src="tutorial.html", style="height:4200px; width:100%", frameborder = "0", scrolling="no")
    helpweb
  })
}

# Run the application
shinyApp(ui = ui, server = server)
