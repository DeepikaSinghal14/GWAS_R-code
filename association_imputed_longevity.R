  ###Aim: Perform association analysis with the imputed data set  
## Deepika Singhal

  library(qqman)
  library(data.table)
  library(tidyverse)
  require(GWASTools)
  library(RColorBrewer)  
  ##### set up env. var ######
  mypath <- if(class(try(rstudioapi::getSourceEditorContext()$path,silent=TRUE))=="try-error") {paste0(dirname(substring(argv[grep("--file=", argv)], 8)),'/')} else paste(head(strsplit(rstudioapi::getSourceEditorContext()$path,'/')[[1]],-1),collapse='/')
      
  ##### Association testing #####
system(paste0("plink --bfile GSA_imp_3103_QCed  --logistic --freq case-control --covar GSAhapmap.eigenvec --covar-number 1-5 --out GSA_imp"))
GSAassocImp <- fread("GSA_imp.assoc.logistic", data.table = F)
GSAassocImp_add <- GSAassocImp[GSAassocImp$TEST == "ADD", ] # 7363285 variants
sigSNPsImp <- GSAassocImp_add %>% arrange(P) %>% filter(P <= 6e-9)
SNPlistrawImp <- fread("GSA_imp_3103_QCed.bim", data.table = F)
sigSNPsImp <- sigSNPsImp %>% left_join(SNPlistrawImp, by=c('BP' = 'V4')) %>% select(-"V1", -"V3",-"V5", -"V6")
SNPshighlightImp <- sigSNPsImp %>% select("SNP")
write.table(SNPshighlightImp,"SNPshighlightImp", col.names = T, row.names = F, sep = " ", quote = F)
      
  

  ### Correcting for multiple testing
0.05/7363285=6e-9 
  
    # calculate manually      
chiSQ <- qchisq(1- GSAassocImp_add$P, 1)
lambda <- round(median(chiSQ, na.rm = T)/qchisq(0.5, 1), digits = 5)
      
  ##Manhattan plot
    GSAassocImp <- GSAassocImp_add[ , c(1,2,3,9)]  
    GSAassocImp <- GSAassocImp[!is.na(GSAassocImp$P), ]
    write.table(GSAassocImp,"GSAassocImp", col.names = T, row.names = F, sep = " ", quote = F)
    manhattan(GSAassocImp, chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: imputed", suggestiveline = F, genomewideline = -log10(6e-9), cex = 0.3, ylim = c(0, 25))
  
 
 
