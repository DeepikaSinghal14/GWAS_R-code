#### Candidate association analsysis ####
# Deepika Singhal
# Creation: 07.04.2020
setwd("/home/deepika/")
library(readODS)
library(tidyverse)
library(data.table)
library(qqman)

### Create SNP lists for subsetting the data set ###
# import data from ods and change names to data format (chr_:_) and write SNP list
candidates <- read_ods(path = "/home/deepika/Documents/fertility_subset.ods")
candidates <- candidates%>%select("chromosome", "bp ")
candidates <- candidates%>%mutate(new_name=paste0("chr",candidates$chromosome,":",candidates$bp))
head(candidates)
candidates <- candidates%>%select("new_name")
write.table(candidates, "candidates_all", row.names = F, col.names = F, sep = " ", quote = F)
candidates <- fread("candidates_all", header = F)

### Check for duplicates ###
dups <- candidates[duplicated(candidates$V1)]
candidates_nodups <- candidates%>%anti_join(dups)
dups_in_dups <- dups[duplicated(dups$V1)]
dups_wo_dups <- dups%>%anti_join(dups_in_dups)
dups_in_dups <- dups_in_dups[duplicated(dups_in_dups$V1)]
candidates_nodups <- candidates_nodups%>%rbind(dups_wo_dups)%>%rbind(dups_in_dups)
dups <- candidates_nodups[duplicated(candidates_nodups$V1)]
candidates_nodups <- candidates_nodups%>%anti_join(dups)
candidates_nodups <- candidates_nodups%>%rbind(dups)
dups <- candidates_nodups[duplicated(candidates_nodups$V1)]
write.table(candidates_nodups, "candidates_all", row.names = F, col.names = F, sep = " ", quote = F)

### Extract candidates ###
system("plink --bfile GSA_imp_3103_QCed --extract candidates_all --make-bed --out subset_fertility")


### Candidate_fertility association testing ###
system(paste0("plink --bfile subset_fertility  --logistic --covar GSAhapmap.eigenvec --covar-number 1-5 --out subset_fertility_menopause"))
GSAassoc <- fread("subset_fertility_menopause.assoc.logistic", data.table = F)
GSAassoc_add <- GSAassoc[ GSAassoc$TEST == "ADD", ] # 666 variants 
sigSNPs <- GSAassoc_add %>% arrange(P) %>% filter(P <= 0.05)      

### Correct for multiple testing (Bonferroni) ###
0.05/666 # 7.507508e-05


### Visualization ###

GSAassoc <- GSAassoc_add[ , c(1,2,3,9)]  
GSAassoc <- GSAassoc[!is.na(GSAassoc$P), ]
write.table(GSAassoc,"GSAassoc", col.names = T, row.names = F, sep = " ", quote = F)
manhattan(GSAassoc, chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: GSA", suggestiveline = F, genomewideline = -log10(7.507508e-05), cex = 0.3, ylim = c(0, 25))
  