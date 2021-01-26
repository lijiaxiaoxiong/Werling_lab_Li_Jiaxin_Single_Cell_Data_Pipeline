# Jan 24th 2021
cell_type <- "AST-PP"

# Get the exprMatrix and meta of this cell type
expr_cell <- get(paste0("expr_", cell_type))
meta_cell <- get(paste0("meta_", cell_type))

# Set the simulation numbers and logs here: 
simulation_number <- 20  # The round number of simulations
log_address <- "Jan_24th_Bootstrap_log.txt"
file_path <- paste0("/mnt/sas0/AD/jli2274/logs/", log_address)
write_lines("", file = file_path) 


# Set some symbols
segmentation <- "*************************************************************"

set.seed(42)
for(i in 1:simulation_number){
  random_list <- sample(control_male_sample_list, 4, replace = F)
  
  sample_list <- c(random_list, control_female_sample_list)
  
  # Prepare the expr and meta
  meta <- dplyr::filter(meta_cell, meta_cell$sample %in% sample_list)
  cell_list <- unlist(meta$cell)
  expr <- dplyr::select(expr_cell, c("gene_ID", cell_list))
  
  # Screen out genes that are under-expressed
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
  zlmCond <- zlm(~sex + cngeneson + region + Capbatch + Seqbatch + `RNA.Integrity.Number` + `RNA.mitochondr..percent` + `RNA.ribosomal.percent` + age, 
                 sca, method = "bayesglm", ebayes = F, silent =T)
  
  
  # Summary the condition
  summaryCond <- summary(zlmCond, doLRT = 'sexM')
  summaryDt <- summaryCond$datatable
  
  
  # Merge the logFC and pvalue together 
  fcH <- merge(summaryDt[contrast=='sexM' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
               summaryDt[contrast=='sexM' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  
  # Calculate the fdr 
  fcH[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  
  # Calculate the 
  fcH[,chi_sq_stat:=qchisq(`Pr(>Chisq)`, df = 1)]
  
  assign(paste0("fcH_",cell_type,"_" ,i), fcH)
  
  # Write information into log
  write_lines(segmentation, file = file_path, append = T)
  write_lines(paste0("Simulation ", i), file = file_path, append = T)
  write_lines("Female sample list:", file = file_path, append = T)
  write_lines(control_female_sample_list, file = file_path, append = T)
  write_lines("Male sample list:", file = file_path, append = T)
  write_lines(random_list, file = file_path, append = T)
  
  cell_number <- length(cell_list)
  write_lines(paste0("Total cell number: ", cell_number), file = file_path, append = T)
  female_cell_number <- dim(meta[meta$sample %in% control_female_sample_list,])[1]
  male_cell_number <- dim(meta[meta$sample %in% random_list, ])[1]
  write_lines(paste0("F:M ",female_cell_number, ":", male_cell_number), file = file_path, append = T)
  print("Simulation %d has been done! Seed: 42", i)
}