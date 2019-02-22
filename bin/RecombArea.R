library(Matrix)
library(blockseg)
library(reshape2)
library(tidyr)
options(warn = -1)
rm(list=ls())

source("/pasteur/homes/aghozlan/entvirus/bin/getBreaks_th.R")
args <- commandArgs(TRUE)
print(args)
print(args[4])
if (length(args)==0) {
  stop("At least three arguments must be supplied (input data file, input simul file and output file).\n", call.=FALSE)
} else if (length(args)>1) {

  ## input data matrix
  file1 <- as.character(args[1])
  ## input simul matrix
  file2 <- as.character(args[2])
  ## output file
  output <- as.character(args[3])
  ## optional threshold for clustering
  if(is.na(args[4]) || args[4]<=0 || args[4]>=100) th_clust = 40 else th_clust =args[4]
}

## Loading
  m2 = as.matrix(read.csv(file = file1, sep="\t", header = FALSE))
  m2_sim = as.matrix(read.csv(file = file2, sep="\t", header = FALSE))

  stopifnot(nrow(m2)==nrow(m2_sim),ncol(m2)==ncol(m2_sim))

## Get the breakpoints
  set.seed(22)
  restab = stab.blockSeg(m2,500,max.break =min(25,trunc(nrow(m2)/2)),sym.break = FALSE,random.break = FALSE,verbose = FALSE,mc.cores=4)

## Get the summary stat for the 2 matrices
  res = getBreaks_th(restab,m2,threshold=th_clust,postprocessing=list(post=TRUE,adjacent=1),plot=FALSE)
  res_sim = getBreaks_th(restab,m2_sim,threshold=th_clust,postprocessing=list(post=TRUE,adjacent=1),plot=FALSE)

## Get the quantiles of the student for each clusters (risk = 1%, multiple correction with bonferonni)
  mat_qt = qt(1-0.005/prod(dim(res$resum)),df=res_sim$emplz%*%t(res_sim$emplw)-1)
  mat_qt[is.na(mat_qt)]=0

## Define the upper and lower bounds of the CI
  simBorneSup = res_sim$mean+mat_qt*sqrt(res_sim$var/(res_sim$emplz%*%t(res_sim$emplw)))
  dataBorneInf =res$mean-mat_qt*sqrt(res$var/(res$emplz%*%t(res$emplw)))

## Check for recomb
  matres = simBorneSup<dataBorneInf

## if NA, ie if only one val, compare value to the max
  matres[is.na(matres)] = (res$mean>max(m2_sim,na.rm = T))[is.na(matres)]

## Add the names to matres
  # columns
  cnames = cumsum(res_sim$emplw)*100
  colnames(matres) = paste(c(0,cnames[-length(cnames)]),cnames,sep="--")
  # rows
  rnames = cumsum(res_sim$emplz)*100
  rownames(matres) = paste(c(0,rnames[-length(rnames)]),rnames,sep="--")

## Keep only the recombined area
  mat_melt = melt(matres,varnames = c("Rows","Cols"),value.name = "Recomb")
  mat_melt = mat_melt[which(mat_melt$Recomb),]

  mat_melt = separate(data=mat_melt,col=Rows,into=c("Row_start","Row_stop"),sep="--")
  mat_melt = separate(data=mat_melt,col=Cols,into=c("Col_start","Col_stop"),sep="--")

## Combine by row
  ind_to_remove_row = c()
  for( i in 1:(nrow(mat_melt)-1))
  {
    if(mat_melt$Row_stop[i]==mat_melt$Row_start[i+1] & mat_melt$Col_start[i]==mat_melt$Col_start[i+1] & mat_melt$Col_stop[i]==mat_melt$Col_stop[i+1] )
    {
      mat_melt$Row_start[i+1] = mat_melt$Row_start[i]
      ind_to_remove_row = c(ind_to_remove_row,i)
    }
  }
  if(length(ind_to_remove_row)>0) out_tmp  = mat_melt[-ind_to_remove_row,]
  out_tmp = out_tmp[order(out_tmp$Row_start),]

## Combine by col
  ind_to_remove_col = c()
  for( i in 1:(nrow(out_tmp)-1))
  {
    if(out_tmp$Col_stop[i]==out_tmp$Col_start[i+1] & out_tmp$Row_start[i]==out_tmp$Row_start[i+1] & out_tmp$Row_stop[i]==out_tmp$Row_stop[i+1] )
    {
      out_tmp$Col_start[i+1] = out_tmp$Col_start[i]
      ind_to_remove_col = c(ind_to_remove_col,i)
    }
  }

  if(length(ind_to_remove_col)>0) out  = out_tmp[-ind_to_remove_col,] else out=out_tmp

## Output
  write.csv(out,file=output,quote=FALSE,col.names = TRUE,row.names=FALSE)
