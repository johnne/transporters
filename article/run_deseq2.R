# Statistical analysis using DESeq2
library("DESeq2")
# Functionality
doy <- function(samples){
  doys <- rep(1,length(samples))
  i = 1
  for (sample in samples){
    x = as.Date(sample,"%y%m%d")
    day = as.numeric(strftime(x, format = "%j"))
    doys[i] = day
    i = i+1
  }
  doys
}

dseq_groups <- function(coldata,countdata,group1,group2){
  coldata = as.data.frame(coldata)
  this_coldata = coldata[coldata$group%in%c(group1,group2),]
  
  this_coldata$group = droplevels(this_coldata$group)
  dseq = DESeqDataSetFromMatrix(countData = countdata[,rownames(this_coldata)],colData = this_coldata, design = ~group)
  dseq_run = DESeq(dseq)
  res = results(dseq_run)
  string=paste(group2,group1,sep=" vs ")
  res$groups = rep(string,nrow(res))
  res$transporter = row.names(res)
  res
}

## Read filtered counts data for abundant transporters
#Because the power of statistical analyses of metagenomes depends heavily on the number of reads mapped to features, 
#only TIGRFAMs with >= 100 average counts across datasets were included here. 
#The representative TIGRFAMs of each transporter was used as a proxy for that transporter and the read counts summed.
mg_counts <- read.table("results/mg/rep_trans_filt.raw_counts.tsv", header=T, sep="\t", row.names=1)
colnames(mg_counts) <- gsub("X","",colnames(mg_counts))

mt_counts <- read.table("results/mt/rep_trans_filt.raw_counts.tsv", header=T, sep="\t", row.names=1)
colnames(mt_counts) <- gsub("X","",colnames(mt_counts))

## Set sample metadata
### Read group information on samples
mg_groups <- read.table("results/mg/samplegroups.tab", header=T, sep="\t", col.names = c("Sample","Group"))
rownames(mg_groups) = mg_groups$Sample
mg_groups[order(mg_groups$Group),]

mt_groups <- read.table("results/mt/samplegroups.tab", header=T, sep="\t", col.names = c("Sample","Group"))
rownames(mt_groups) = mt_groups$Sample
mt_groups[order(mt_groups$Group),]

#Use only samples that are clustered into groups.
mg_counts <- mg_counts[as.character(mg_groups$Sample)]
mt_counts <- mt_counts[as.character(mt_groups$Sample)]

### Create dataframe with sample groups
mg_coldata <- data.frame(row.names = colnames(mg_counts), sample=colnames(mg_counts), group=mg_groups[colnames(mg_counts),"Group"],day=doy(colnames(mg_counts)))
mt_coldata <- data.frame(row.names = colnames(mt_counts), sample=colnames(mt_counts), group=mt_groups[colnames(mt_counts),"Group"])

## Analysis of the MG dataset
groups <- levels(mg_groups$Group)
done <- c()
df <- data.frame(baseMean=NA,log2FoldChange=NA,lfcSE=NA,stat=NA,pvalue=NA,padj=NA,groups="",transporter="")
for (group1 in groups){
  for (group2 in groups[!(groups%in%c(group1))]) {
    sort_groups = sort(c(group2,group1))
    x = paste(sort_groups[1],sort_groups[2],sep="X")
    if (!x%in%done){
      df_dseq = as.data.frame(dseq_groups(coldata = mg_coldata, 
                                          countdata = mg_counts, group1=group1,group2=group2))
      df <- rbind(df,df_dseq)
      done <- c(done,x)
    }
  }
}
mg_df = df[2:nrow(df),]
# Adjust p-values
mg_df$p_total_adj = p.adjust(mg_df$pvalue,method="fdr")
write.table(x=mg_df,"results/mg/deseq2.tab",sep="\t",row.names = FALSE, quote=FALSE)

## Analysis of the MT dataset
groups <- levels(mt_groups$Group)
done <- c()
df <- data.frame(baseMean=NA,log2FoldChange=NA,lfcSE=NA,stat=NA,pvalue=NA,padj=NA,groups="",transporter="")
for (group1 in sort(groups)){
  for (group2 in groups[!(groups%in%c(group1))]) {
    sort_groups = sort(c(group2,group1))
    x = paste(sort_groups[1],sort_groups[2],sep="X")
    if (!x%in%done){
      df_dseq = as.data.frame(dseq_groups(coldata = mt_coldata, 
                                          countdata = mt_counts, group1=group1,group2=group2))
      df <- rbind(df,df_dseq)
      done <- c(done,x)
    }
  }
}
mt_df = df[2:nrow(df),]
mt_df$p_total_adj = p.adjust(mt_df$pvalue,method="fdr")
write.table(x=mt_df,"results/mt/deseq2.tab",sep="\t",row.names = FALSE, quote=FALSE)
