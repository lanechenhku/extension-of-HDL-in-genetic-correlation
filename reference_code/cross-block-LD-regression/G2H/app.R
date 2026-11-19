#----
effect_size_plot<-function(sources='GWAS ATLAS',
                           SNPid1=NA, 
                           SNPid2=NA,
                           SNPid3=NA,
                           SNPid4=NA,
                           SNPid5=NA,
                           path = "C:/Users/22367/Downloads",
                           SNPpair=SNPpair,
                           chosen='-log Pval',
                           logbase=10){
  #py_module_available("pandas")#检查“pandas”是否安装
  #path <- "C:/Users/22367/Downloads"
  if(sources=='GWAS ATLAS'){
    file_default<-'Genome wide association study ATLAS'}
  
  fileNames <- dir(path) 
  fileN=fileNames[grep(file_default,fileNames)]
  file.remove(paste0(path,'/',fileN))
  
  #input----
  # SNPid1=NA
  # SNPid2=NA
  # SNPid3=NA
  # SNPid4=NA
  # SNPid5=NA
  
  
  #logbase=10
  SNPlist=na.omit(c(SNPid1,SNPid2,SNPid3,SNPid4,SNPid5))
  SNPno=length(SNPlist)
  #xpath='/html/body/div[1]/div[2]/div/div[4]/div[2]/div/div/div[1]/a[2]'
  if(sources=='GWAS ATLAS'){
    URL=c("https://atlas.ctglab.nl/PheWAS")
  }
  write.table(c(SNPpair,URL),"C:/Users/22367/Desktop/thesis/cross-block-LD-regression/old G2H/r.txt",row.names = F,col.names = F,quote = F)
  
  #get data----
  
  source_python("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/getdata.py")
  # py_run_file("C:/Users/22367/Desktop/thesis/cross-block-LD-regression/old G2H/getdata.py")
  
  #read file----
  path <- "C:/Users/22367/Downloads"
  #aaa=paste0('Genome wide association study ATLAS (',c(1:(SNPno-1)),').csv')
  fileNames <- c('Genome wide association study ATLAS.csv','Genome wide association study ATLAS (1).csv')
  filePath <- sapply(fileNames, function(x){ 
    paste(path,x,sep='/')})   
  data <- lapply(filePath, function(x){
    read.csv(x)}) 
  
  #----
  
  # id1=which(SNPlist==SNPpair[1])
  # id2=which(SNPlist==SNPpair[2])
  # df1=data[[id1]]
  # df2=data[[id2]]
  df1=data[[1]]
  df2=data[[2]]
  df=left_join(df1,df2,by='atlas.ID')
  df=na.omit(df)
  variable=c("atlas.ID",  "PMID.x" ,"Domain.x" , "Trait.x", "P.value.x","P.value.y")
  dd=df[,variable]
  # dd$z.x=as.numeric(dd$beta.x)/as.numeric(dd$sd.x)
  # dd$z.y=as.numeric(dd$beta.y)/as.numeric(dd$sd.y)
  if(chosen=='-log Pval'){
    dd$z.x=-log(as.numeric(dd$P.value.x),logbase)
    dd$z.y=-log(as.numeric(dd$P.value.y),logbase)
  }
  dd$dataset=as.character(dd$PMID.x)
  dd$group=dd$Domain.x
  dd$dataset[grep('BioRxiv',dd$PMID.x)]='BioRxiv'
  # dd$dataset[grep('GWAS_CKDGenConsortium',dd$dataset.x)]='GWAS_CKDGenConsortium'
  # dd$dataset[grep('GWAS_GEFOS',dd$dataset.x)]='GWAS_GEFOS'
  # dd$dataset[grep('GWAS_GODARTS',dd$dataset.x)]='GWAS_GODARTS'
  # dd$dataset[grep('GWAS_HOORN',dd$dataset.x)]='GWAS_HOORN'
  # dd$dataset[grep('TOPMED',dd$dataset.x)]='TOPMED'
  p1<-ggplot(data=dd,aes(x=z.x,y=z.y,col=dataset))+geom_point()+xlab(SNPpair[1])+ggtitle(paste0(chosen,' scatterplot'))+
    geom_line(stat='smooth',method = 'lm', formula = y ~ x+1,se=F,alpha=0.5,size=0.5) + #绘制回归直线
    stat_cor(data = dd, method = "pearson")
  
  p2<-ggplot(data=dd,aes(x=z.x,y=z.y,col=group))+geom_point()+xlab(SNPpair[1])+ylab(SNPpair[2])+ggtitle(paste0(chosen,' scatterplot'))+
    geom_line(stat='smooth',method = 'lm', formula = y ~ x+1,se=F,alpha=0.5,size=0.5) + #绘制回归直线
    stat_cor(data = dd, method = "pearson")
  p<-list(p1,p2)
  return(p)
}
#----

library(reticulate)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
library(shiny)
library(shinythemes)
library(shinyjs)
library(knitr)
library(reshape2)
library(plotly)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(DT)
library(shinydashboard)
library(plotly)
library(RColorBrewer)
library(ggpmisc)
library(Rmisc)

library(reticulate)
library(tidyverse)
library(ggpubr)
library(ggpmisc)
SNPid1=c('rs78378222')
SNPid2=c('rs3730556')
SNPpair=c(SNPid1,SNPid2)

choice=c('-log Pval','Beta','Z Score')
database=c('GWAS ATLAS','GWAS Catalog','MSKPP')
#----

database=c('GWAS ATLAS','msk.hugeamp.org','GWAS catalog')
choicee=c('-log Pval','Beta','Z score')
#main structure of webpage----
ui <- dashboardPage(
    dashboardHeader(title = "Bi-variate Summary statistics Profile of SNP Pair Across GWASs"),
    
    dashboardSidebar(
        selectizeInput(inputId ='sources', 
                       label = "Data sources", 
                       choices = as.factor(database),
                       options = list(
                           placeholder = 'Choose a database'
                       ),
                       selected='GWAS ATLAS'
                       
        ),
        selectizeInput(inputId ='chosen', 
                       label = "type of summary statistics", 
                       choices = as.factor(choice),
                       options = list(
                           placeholder = 'Please choose one'
                       ),
                       selected='-log Pval'
                       
        ),
        sliderInput("logbase",
                    "Base of log (if -logP is chosen)",
                    min = 2,
                    max = 10,
                    step = 0.01,
                    value = 2.72),
        textInput(inputId='SNPid1', 'SNP1: rsID', value = 'rs78378222'),
        textInput(inputId='SNPid2', 'SNP2: rsID', value = 'rs3730556'),
        textInput(inputId='SNPid3', 'SNP3: rsID', value = 'rs9894946'),
        textInput(inputId='SNPid4', 'SNP4: rsID', value = 'rs78378222'),
        textInput(inputId='SNPid5', 'SNP5: rsID', value = 'rs78378222'),
        tags$p('At least input two SNPs')),
#主体函数

     dashboardBody(
       
        # fluidRow(
        #     valueBoxOutput("traits"),
        #     valueBoxOutput("high_risk"),
        #     valueBoxOutput("low_risk")
        # ),
        # fluidRow(
        #     box(title = "Relative Genetic Risk Profile [RGR = (your GR/population average GR-1)%]", 
        #         subtitle='aa',
        #         solidHeader = T,
        #         width = 15, collapsible = T,
        #         plotlyOutput("grisk_plot"))
        # ), #row
        # fluidRow(
        #     box(title = "Estimated heritability, onset age and population prevalence for the disease category", solidHeader = T,
        #         width = 15, collapsible = T,
        #         plotlyOutput("grisk_plot3"))
        # ), #row
        fluidRow(
            box(title = "Stratify by disease groups", solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("p1"))
        ), #row
        fluidRow(
            box(title = "Stratify by datasets (publications）", solidHeader = T,
                width = 15, collapsible = T,
                plotlyOutput("p2"))
        ) #row
    ) # body
    
)


server <- function(input, output) {
  sources=reactive(input$sources)
  SNPid1=reactive(input$SNPid1)
  SNPid2=reactive(input$SNPid2)
  SNPid3=reactive(input$SNPid3)
  SNPid4=reactive(input$SNPid4)
  SNPid5=reactive(input$SNPid5)
 # path = reactive(input$path)
  SNPpair=reactive(c(input$SNPid1,input$SNPid2))
  chosen=reactive(input$chosen)
  logbase=reactive(input$logbase)

  
  output$p1 <- renderPlotly({
    sources=sources()
    SNPid1=SNPid1()
    SNPid2=SNPid2()
    SNPid3=SNPid3()
    SNPid4=SNPid4()
    SNPid5=SNPid5()
   # path = path()
    SNPpair=SNPpair()
    chosen=chosen()
    logbase=logbase()
    p<-effect_size_plot(
    sources=sources,
    SNPid1=SNPid1,
    SNPid2=SNPid2,
    SNPid3=SNPid3,
    SNPid4=SNPid4,
    SNPid5=SNPid5,
    #path = path,
    SNPpair=SNPpair,
    chosen=chosen,
    logbase=logbase)
   
        ggplotly(p[[1]])
    })
    
    output$p2 <- renderPlotly({
        ggplotly(p[[2]])
    })
    
    
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)
