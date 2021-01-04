###Aim: Perform association analysis with unimputed data set  
#### set up env. var ######
mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
setwd(mypath)

##### Association testing #####
system(paste0("plink --bfile GSAdata_QCed_3003_assoc  --logistic --freq case-control --covar GSAhapmap.eigenvec --covar-number 1-5 --out GSAunimp"))
GSAassoc <- fread("GSAunimp.assoc.logistic", data.table = F)
GSAassoc_add <- GSAassoc[ GSAassoc$TEST == "ADD", ] 
sigSNPs <- GSAassoc_add %>% arrange(P) %>% filter(P <= 1e-05)
SNPlistraw <- fread("GSAdata_QCed_3003_assoc.bim", data.table = F)
sigSNPs <- sigSNPs %>% left_join(SNPlistraw, by=c('BP' = 'V4')) %>% select(-"V1", -"V3",-"V5", -"V6") 
SNPshighlight <- sigSNPs %>% select("SNP")
write.table(SNPshighlight,"SNPshighlight", col.names = T, row.names = F, sep = " ", quote = F)


##Manhattan plot

GSAassoc <- GSAassoc_add[ , c(1,2,3,9)]  
GSAassoc <- GSAassoc[!is.na(GSAassoc$P), ]
write.table(GSAassoc,"GSAassoc", col.names = T, row.names = F, sep = " ", quote = F)
manhattan(GSAassoc, chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: GSA", col=brewer.pal(3, "Set2"), cex = 0.3, ylim = c(0, 25))


