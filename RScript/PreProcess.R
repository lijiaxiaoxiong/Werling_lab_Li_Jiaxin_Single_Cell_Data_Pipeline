# Decide the cell type
cell_type <- "L4"
# read cell type meta
meta <- NULL
expr <- NULL
read_data <-function(cell_type){
  dataroot <- "/Users/apple/Box/Li_Jiaxin/SingleCellASDBrain/data/Cell_type/"
  meta_path <- paste0(dataroot, cell_type, "_meta.csv")
  expr_path <- paste0(dataroot, cell_type, "_expr.csv")
  meta <<- read.csv(meta_path)
  print(cell_type)
  print(" meta file load!")
  expr <<- read.csv(expr_path)
  print(" expr file load!")
}
read_data(cell_type)



# Delete the index added
delete_x <- function(df){
  df <- dplyr::select(df, -c("X"))
  return(df)
}

#meta <- delete_x(meta)
#expr <- delete_x(expr)


meta <- `meta_AST-PP_ASD`
expr <- `expr_AST-PP_ASD`

# Filter the genes which is not expressed in over 90% cells
filter_genes <- function(exprMatrix, percentage){
  only_num <- dplyr::select(exprMatrix, -c("gene_ID"))
  exprMatrix$filter_cri <- rowMeans(only_num > 0) >= percentage
  filtered_expr <- filter(exprMatrix, filter_cri == T)
  filtered_expr <- dplyr::select(filtered_expr, -("filter_cri"))
  return(filtered_expr)
}

filted_expr <- filter_genes(expr, 0.1)
gene_list <- filted_expr$gene_ID
gene_listdf <- data.frame(primerid = gene_list)

fdata <- gene_listdf
cdata <- meta
names(cdata)[1] <- "wellKey"
exprMatrix <- as.matrix(dplyr::select(filted_expr, -c("gene_ID")))


# Make the sca object
sca <- MAST::FromMatrix(exprMatrix, cData = cdata, fData = fdata)


# Calculate the cngeneson
cdr2 <- colSums(assay(sca) > 0)
colData(sca)$cngeneson <- scale(cdr2)


# Building the model
zlmCond <- zlm(~sex + cngeneson + diagnosis + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
               sca, method = "bayesglm", ebayes = F, silent =T)

# Building the model
zlmCond <- zlm(~sex + cngeneson + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
               sca, method = "bayesglm", ebayes = F, silent =T)


# Summary the condition
summaryCond <- summary(zlmCond, doLRT = 'sexM')
summaryDt <- summaryCond$datatable


# Merge the logFC and pvalue together 
fcH <- merge(summaryDt[contrast=='sexM' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
              summaryDt[contrast=='sexM' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients

`fcH_AST-PP_ASD` <- fcH

# Calculate the fdt filter out the DEGs
fcH[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
DEGs <- merge(fcH[fdr<.05 & abs(coef)> log2(1.1)], as.data.table(mcols(sca)), by='primerid')



# Set the plot path
plotpath <- paste0("/Users/apple/Box/Li_Jiaxin/SingleCellASDBrain/results/", cell_type, "/Plot/")


# Volcano Plot
png(filename = paste0(plotpath, "volcano.png"))
ggplot(DEGs, aes(coef, -log10(fdr))) + geom_point() + ggtitle("Volcano After Filter") + xlab("log2 FC")  + ylab("-log10 fdr") + theme_light()
dev.off()


make_gene_profile <- function(gene_list, DEGs){
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  bm_result <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"), 
                     filters = "ensembl_gene_id",values = gene_list, mart = ensembl)
  names(bm_result)[1] <- "primerid"
  gene_profile <- merge(DEGs, bm_result, all = T, by <- "primerid")
  return(gene_profile)
}

DEG_profile <- make_gene_profile(DEGs$primerid, DEGs)

file_path <-  paste0("/Users/apple/Box/Li_Jiaxin/SingleCellASDBrain/results/", cell_type, "/")
write.csv(DEG_profile, file = paste0(file_path, "DEGs.csv"))

DEG_profile_bad <- DEG_profile[DEG_profile$coef * DEG_profile$ASD_logFC < 0 & DEG_profile$coef * DEG_profile$Control_logFC < 0,]
