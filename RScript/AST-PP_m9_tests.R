# Jan 19th
# Author: Jiaxin Li
fcH_AST_PP_diagnosis_m9 <-`fcH_AST-PP_diagnosis`[`fcH_AST-PP_diagnosis`$primerid %in% Pari_m9 ]

sig_diagnosis_m9 <- fcH_AST_PP_diagnosis_m9[fcH_AST_PP_diagnosis_m9$fdr < 0.05,]

ggplot(fcH_AST_PP_diagnosis_m9, aes(y = coef)) + geom_boxplot()

ggplot(sig_diagnosis_m9, aes(y = coef)) + geom_boxplot()


# t-test of m9 ASD and control

qq1 <- qqnorm(fcH_AST_PP_m9$coef.x, main = "ASD logFC Q-Q Plot") #+ qqline(fcH_AST_PP_m9$coef.x, col = 2, lwd = 2)
qq2 <- qqnorm(fcH_AST_PP_m9$coef.y, main = "Control logFC Q-Q Plot")


t.test(fcH_AST_PP_m9$coef.x, fcH_AST_PP_m9$coef.y, paired = T)

wilcox.test(fcH_AST_PP_m9$coef.x, fcH_AST_PP_m9$coef.y, paired = T)
s
