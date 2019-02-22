library(plyr)
library(data.table)
# input=function(){
#   x=commandArgs(trail=T)
#   if(length(x) != 3){
#     #id% cov%
#     print("Usage : ./virus_overtime.R result_summary.tsv  sample_target_v5 groups")
#     q("yes",status=1)
#   }
#   return(x)
# }
# x = input()
# 
# resume = read.table(x[1], sep="\t")
# target = read.table(x[2], sep="\t")
# groups = unlist(strsplit(scan(x[3], what="character", sep=NULL)[3],","))


resume = read.table("result_summary.tsv", sep="\t", header=T)
#groups = unlist(strsplit(scan("group_aa", what="character", sep=NULL)[3],","))
target = read.table("sample_target_all_v5.csv", sep="\t",header=T)
groups = strsplit(as.character(read.delim("groups.txt", sep="\t", header=F)[,3]), ",")

# take each group
group_num = 1
for (vp1 in groups){
  sample_names = tstrsplit(vp1, "|", fixed=T)[[1]]
  # associate to each vp1 a date
  collect = sapply(sample_names, function(x) target[which(x == target$sample),]$COLLECTION_DATE)
  
  # associate to each vp1 an annotation
  serotype = resume[which(resume$VP1_sequences %in%vp1),]$Serotype_VP1
  specie = resume[which(resume$VP1_sequences %in%vp1),]$Specie_VP1
  
  # group all
  rah = data.frame(vp1,as.Date(collect), serotype, specie)
  ser = unique(rah$serotype)
  spe = unique(rah$specie)
  # Count per month
  grouped = rah %>% group_by(month=floor_date(as.Date(collect), "month")) %>% tally()
  svg(paste(paste(spe,ser,"group",group_num,sep="_"), ".svg",sep=""),width=21,height=7)
  plot(grouped,xaxt="n", type="l", main=paste(spe,ser))
  points(grouped)
  axis(1, grouped$month, format(grouped$month, "%b-%y"), cex.axis = .7)
  group_num = group_num +1
  dev.off()
}

