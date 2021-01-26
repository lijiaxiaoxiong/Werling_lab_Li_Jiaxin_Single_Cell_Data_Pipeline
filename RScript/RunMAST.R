# Jan 18th 2021
# Author: Jiaxin Li
# A Function to do MAST zlm analyze from expr and meta

filter_genes <- function(exprMatrix, percentage){
  only_num <- dplyr::select(exprMatrix, -c("gene_ID"))
  exprMatrix$filter_cri <- rowMeans(only_num > 0) >= percentage
  filtered_expr <- filter(exprMatrix, filter_cri == T)
  filtered_expr <- dplyr::select(filtered_expr, -("filter_cri"))
  return(filtered_expr)
}

zlm_MAST <- function(expr_with_gene_ID, meta, if_diagnosis_in_model = T){
  expr_f <- filter_genes(expr_with_gene_ID, 0.1)
  gene_list_f <- expr_f$gene_ID
  gene_listdf_f <- data.frame(primerid = gene_list_f)
  
  fdata <- gene_listdf_f
  cdata <- meta
  names(cdata)[1] <- "wellKey"
  exprMatrix <- as.matrix(dplyr::select(expr_f, -c("gene_ID")))
  
  
  # Make the sca object
  sca <- MAST::FromMatrix(exprMatrix, cData = cdata, fData = fdata)
  
  
  # Calculate the cngeneson
  cdr2 <- colSums(assay(sca) > 0)
  colData(sca)$cngeneson <- scale(cdr2)
  
  
  # Building the model
  if(if_diagnosis_in_model){
    zlmCond <- zlm(~sex + cngeneson + diagnosis + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
                   sca, method = "bayesglm", ebayes = F, silent =T)
  }else{
    zlmCond <- zlm(~sex + cngeneson + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
                   sca, method = "bayesglm", ebayes = F, silent =T)
  }
  
  
  
  # Summary the condition
  summaryCond <- summary(zlmCond, doLRT = 'sexM')
  summaryDt <- summaryCond$datatable
  
}