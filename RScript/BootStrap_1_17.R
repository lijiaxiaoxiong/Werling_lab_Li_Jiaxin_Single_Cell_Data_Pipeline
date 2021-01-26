# Jan 17th
# Author: Jiaxin Li
# --------------------------------------------------------
# 1. Sample Stat
# Sample stats
sample_list <- unique(meta_whole$sample)
sample_stat <- data.frame(sample = sample_list, cell_number = 0)# The cell_number column need to be reassigned as numbers, so they should to be numbers at the beginning.sy
for(s in sample_list){
  meta_sample <- dplyr::filter(meta_whole, meta_whole$sample == s)
  print(dim(meta_sample))
  sample_stat$cell_number[sample_stat$sample == s] <- dim(meta_sample)[1] # 
}
  
# What happened to this? why I need to initiate the dataframe firstly?

meta_control_female <- filter(meta_whole, meta_whole$diagnosis == "Control" & meta_whole$sex == "F")

control_female_sample_list <- unique(meta_control_female$sample)

control_female_sample_stat <- data.frame(sample = control_female_sample_list, cell_number = 0)# The cell_number column need to be reassigned as numbers, so they should to be numbers at the beginning.sy
for(s in sample_list){
  meta_sample <- dplyr::filter(meta_control_female, meta_control_female$sample == s)
  print(dim(meta_sample))
  control_female_sample_stat$cell_number[control_female_sample_stat$sample == s] <- dim(meta_sample)[1] # 
}


meta_control_male <- filter(meta_whole, meta_whole$diagnosis == "Control" & meta_whole$sex == "M")
control_male_sample_list <- unique(meta_control_male$sample)

list_log <- c()
# Randomly select male sample set
set.seed(42)
`meta_AST-PP` <- filter(meta_whole, meta_whole$cluster == "AST-PP")

simulation_log <- data.frame(ID = c(1:10), Male_samples = rep("a", 10), Male_cell_number = c(1:10), Female_cell_number = c(1:10))


random_list <- sample(control_male_sample_list, 4, replace = FALSE)
list_log <- c(list_log, random_list)
meta_male <- filter(`meta_AST-PP`, `meta_AST-PP`$sample %in% random_list)
  #
meta_female <- filter(`meta_AST-PP`, `meta_AST-PP`$sample %in% control_female_sample_list)
  
meta_random <- rbind(meta_male, meta_female)
random_cell_list <- meta_random$cell
expr_random <- dplyr::select(expr_whole, c("gene_ID",random_cell_list))

filtered_expr_random <- filter_genes(expr_random, 0.1)

gene_list_random <- filtered_expr_random$gene_ID
gene_listdf_random <- data.frame(primerid = gene_list_random)

fdata <- gene_listdf_random
cdata <- meta_random
names(cdata)[1] <- "wellKey"
exprMatrix_random <- as.matrix(dplyr::select(filtered_expr_random, -c("gene_ID")))


# Make the sca object
sca <- MAST::FromMatrix(exprMatrix_random, cData = cdata, fData = fdata)


# Calculate the cngeneson
cdr2 <- colSums(assay(sca) > 0)
colData(sca)$cngeneson <- scale(cdr2)


# Building the model
zlmCond <- zlm(~sex + cngeneson  + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
               sca, method = "bayesglm", ebayes = F, silent =T)


# Summary the condition
summaryCond <- summary(zlmCond, doLRT = 'sexM')
summaryDt <- summaryCond$datatable

fcH <- merge(summaryDt[contrast=='sexM' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
             summaryDt[contrast=='sexM' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients


# Calculate the fdt filter out the DEGs
fcH[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
DEGs <- merge(fcH[fdr<.05 & abs(coef)> log2(1.1)], as.data.table(mcols(sca)), by='primerid')

DEG_profile <- make_gene_profile(DEGs$primerid, DEGs)

# random_generator <- function(whole_meta, )

# Jan 24th 2021



