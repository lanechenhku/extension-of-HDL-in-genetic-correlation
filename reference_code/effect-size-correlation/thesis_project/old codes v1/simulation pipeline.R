#packages used----
if(!("devtools" %in% rownames(installed.packages()))) install.packages("devtools"); .rootdir <- "/home2/groups/pcsham/projects/GWASALTAS/Tim"; source("/home2/groups/pcsham/projects/GWASALTAS/Tim/Tmisc/inst/Startup.R")
library(devtools) # Tools to Make Developing R Packages Easier
#library(annotater) # Annotate Package Load Calls
library(RcppArmadillo) # not installed on this machine
library(data.table) # Extension of `data.frame` # Extension of `data.frame`
library(Matrix) # Sparse and Dense Matrix Classes and Methods
library(tidyverse) # Easily Install and Load the 'Tidyverse'
#library(lsgl) # not installed on this machine
#library(lmridge) # not installed on this machine
library(data.table) # Extension of `data.frame` # Extension of `data.frame`
#library(glmnet) # Lasso and Elastic-Net Regularized Generalized Linear Models
#library(KRIS) # not installed on this machine
#library(monomvn) # Estimation for MVN and Student-t Data with Monotone Missingness
#library(DMwR) # Functions and data for "Data Mining with R"


#load functions used by Tim
Tmisc()
Tim.load("Rplink",lib.dir="/home/unfated/Tim")
Tim.load("UKBBscripts", "/psychipc01/disk3/data/UKBB")
Tim.load("lassosum","/home/unfated/Tim")
Tim.load("crosspred", "/home/unfated/Tim")



# get genotype from chr i bfiles and transform into matrix----

input=c('/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train')
bed=paste0(input,i,'.bed')
bim=paste0(input,i,'.bim')
fam=paste0(input,i,'.fam')
block=read.bed(bed,bim,fam)
genodata=block$snp


#get chr 22 genotype directly from txt file ----


genodata=fread('/home/unfated/simulation/chr22_geno_20k_complete.txt') %>%  as.data.frame(data4)


#extract genotype subsets from genotype matrix----

#n.row, n.col: vector for subset to input into LD pruning
n.row=sample(1:nrow(genodata),100000)
n.col=c(1:10000)
#n.ind, n.snp:vector for subset to output
n.ind=c(1:10000)
n.snp=c(1:1000)

extract_geno<-function(genodata,n.row,n.col,threshold,n.ind,n.snp){
  fulldata=genodata
  #ncol(fulldata)
  #nrow(fulldata)

  useddata = fulldata[,which(colSums(abs(fulldata)) !=0)][n.row,n.col][,which(colSums(abs(useddata)) !=0)]
  sum(is.na(useddata))
  #ncol(useddata)
  #nrow(useddata)
  #DO LD pruning so that all LD <=threshold
  useddata1=useddata
  j=1
  while(j<ncol(useddata1)){
    LD=cor(useddata1[,j],useddata1[,-j])
    if(length(which(LD> threshold))>0){
      useddata1=useddata1[,-j]}else{j=j+1}
  }
  
  #ncol(useddata1)
  #nrow(useddata1)
  #sum(is.na(useddata1))
  #remove snps with MAF=0 (genotype constant)
  var.data=apply(useddata1,2,var) 
  useddata1=useddata1[,-which(var.data==0)]
  output.geno=useddata1[n.ind,n.snp]
  return(output.geno)}


#generate phenotype using matrix product----


generate_pheno<-function(output.geno,prop,h2){
  useddata=output.geno
  p=ncol(useddata) 
  #prop: proportion of causal SNPs
  n.causal = ceiling(p*prop)
  beta <- rep(0, p)
  #point normal effect size
  beta[sample(p, n.causal)] <- rnorm(n.causal)
  
  #matrix product
  Xb <- as.matrix(useddata)%*%beta
  #large scale: fast way in paralel using Lassosum

  #heritability
  herit <- h2
  error <- rnorm(length(Xb), sd=sqrt(var(Xb)/herit*(1-herit)))
  y <- Xb + error
  length(y)  
  return(list(y,beta))}

#simulate phenotype in large scale  without using matrix product explicitly----
#load bfile: UKBB whole chromosome 22 genotype

bfile<- "/psychipc01/disk3/data/UKBB/data/genotype/whiteBritish1_train/whiteBritish1_train22"

#check sample size
N<-nrow.bfile(bfile)
#check snp number
p<-ncol.bfile(bfile)


filename <- Rfilename("sim.10000.1000", seed = "prompt")
1
cl <- startParallel(8)
fixed.mem <- 40000

n.causal <-  ceiling(p*prop)

#keep n individuals from all samples
tokeep <- as.logical(sample(N,n))
#use all p SNPs
beta <- rep(0, p)
beta[sample(p, n.causal)] <- rnorm(n.causal)
# use pgs function in Lassosum, after specifying the beta vector
Xb <- pgs(bfile=bfile, weights=beta, cluster=cl, trace=1, 
          keep=tokeep)
herit <- h2
error <- rnorm(n, sd=sqrt(var(Xb)/herit*(1-herit)))
y <- Xb + error
Y <- rep(NA, N)
Y[tokeep] <- y

#get summary statistics using plink: bfile and Y----
sum <- plink(bfile=bfile, cmd="--linear", plink2=TRUE, 
             pheno=Y, silent=F, keep=tokeep, threads = length(cl), 
             mem=fixed.mem)

#or specify files and with covariates: pheno.file, covar.file, pheno, covar
sum <- plink(bfile=bfile, cmd=paste0("--glm"," --pheno ",pheno.file," --pheno-name ",pheno," --covar ",covar.file," --covar-name ",covar), plink2=T, silent=F)

summ <- read.table.plink(sum, ".PHENO1.glm.linear", header=T)
#estimated beta in GWAS
est.beta <- summ$BETA
est.beta[is.na(est.beta)] <- 0
#pval in GWAS
pval <- summ$P
pval[is.na(pval)] <- 1