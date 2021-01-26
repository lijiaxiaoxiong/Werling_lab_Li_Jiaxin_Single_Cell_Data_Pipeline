# Jan 18th 
# Author: Jiaxin Li
# To obseve the M9 patterns in AST-PP
`cell_list_AST-PP` <- unlist(`meta_AST-PP`$cell)
`expr_AST-PP` <- dplyr::select(expr_whole, c("gene_ID") | all_of(`cell_list_AST-PP`))

`meta_AST-PP_control` <- dplyr::filter(`meta_AST-PP`, `meta_AST-PP`$diagnosis == "Control")
`cell_list_AST-PP_control` <- unlist(`meta_AST-PP_control`$cell)
`expr_AST-PP_control` <- dplyr::select(`expr_AST-PP`, c("gene_ID") | all_of(`cell_list_AST-PP_control`))

`meta_AST-PP_ASD` <- dplyr::filter(`meta_AST-PP`, `meta_AST-PP`$diagnosis == "ASD")
`cell_list_AST-PP_ASD` <- unlist(`meta_AST-PP_ASD`$cell)
`expr_AST-PP_ASD` <- dplyr::select(`expr_AST-PP`, c("gene_ID") | all_of(`cell_list_AST-PP_ASD`))


fcH_AST_PP_control_m9 <- dplyr::filter(`fcH_AST-PP_control`, `fcH_AST-PP_control`$primerid %in% Pari_m9)
fcH_AST_PP_ASD_m9 <- dplyr::filter(`fcH_AST-PP_ASD`, `fcH_AST-PP_ASD`$primerid %in% Pari_m9)

fcH_AST_PP_m9 <- merge(fcH_AST_PP_ASD_m9, fcH_AST_PP_control_m9, by = "primerid")
fcH_AST_PP_m9_long <- rbind(fcH_AST_PP_ASD_m9, fcH_AST_PP_control_m9)

ggplot(fcH_AST_PP_m9_long, aes(x = diagnosis, y = coef, color = diagnosis)) + geom_boxplot(colour = "black") + geom_jitter(alpha = 0.5) + geom_line(aes(group = primerid), colour = "gray", alpha = 0.3)+ theme_light() + labs(title = "logFC Comparision between ASD and Control (AST-PP; M9)", y = "logFC")



fcH_AST_PP_ASD_m9[,fdr:=p.adjust(pvalue, 'fdr')]
fcH_AST_PP_control_m9[,fdr:=p.adjust(pvalue, 'fdr')]

sig_control <- fcH_AST_PP_control_m9[fcH_AST_PP_control_m9$fdr < 0.05]
sig_ASD <- fcH_AST_PP_ASD_m9[fcH_AST_PP_ASD_m9$fdr < 0.05]

sig_long <- rbind(sig_ASD, sig_control)
ggplot(sig_long, aes(x = diagnosis, y = coef, color = diagnosis)) + geom_violin() + geom_jitter(alpha = 0.5) + geom_line(aes(group = primerid), colour = "gray", alpha = 0.3)+ theme_light() + labs(title = "logFC Comparision between ASD and Control (AST-PP; M9; Significant)", y = "logFC")



ggplot(fcH_AST_PP_control_m9, aes(y = coef)) + geom_boxplot() + geom_jitter()
ggplot(fcH_AST_PP_ASD_m9, aes(y = coef)) + geom_boxplot()
ggplot(sig_control, aes(y = coef)) + geom_boxplot()
ggplot(sig_ASD, aes(y = coef)) + geom_boxplot()

colnames(`fcH_AST_PP_control_m9`)[3]<- "coef"
colnames(`fcH_AST_PP_ASD_m9`)[3]<- "coef"

colnames(`fcH_AST_PP_control_m9`)[2]<- "pvalue"
colnames(`fcH_AST_PP_ASD_m9`)[2]<- "pvalue"

fcH_AST_PP_m9 <- merge(fcH_AST_PP_ASD_m9, `fcH_AST-PP_control`, by = "primerid")

fcH_AST_PP_m9$diff <- fcH_AST_PP_m9$coef_control - fcH_AST_PP_m9$coef_ASD

sig_short <- fcH_AST_PP_m9[fcH_AST_PP_m9$pvalue_ASD < 0.05 & fcH_AST_PP_m9$`Pr(>Chisq)` < 0.05,]

dim(fcH_AST_PP_m9[fcH_AST_PP_m9$diff > 0,])

fcH_AST_PP_ASD_m9$diagnosis <- "ASD"
fcH_AST_PP_control_m9$diagnosis <- "Control"
fcH_AST_PP_m9_long <- rbind(fcH_AST_PP_ASD_m9, fcH_AST_PP_control_m9)

sig_long <- fcH_AST_PP_m9_long[fcH_AST_PP_m9_long$pvalue < 0.05]

ggplot(fcH_AST_PP_m9_long, aes(x = diagnosis, y = coef)) + geom_boxplot() + geom_jitter() + geom_line(aes(group = primerid))


ggplot(sig_long, aes(x = diagnosis, y = coef, color = diagnosis)) + geom_violin() + geom_jitter(alpha = 0.5) + geom_line(aes(group = primerid, alpha = 0.1), colour = "gray") + theme_light()


