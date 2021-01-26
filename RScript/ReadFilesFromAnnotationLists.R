# Read files from annotation lists:
#cell_type <- "Neu-NRGN-II"
#celltype <-"Neu-NRGN-II"
library(readxl)
# 1. Read Werling 2016
anno <- "/"

# Backgrounds
adultback1 <- read.table(paste0(anno, "werling2016_adult1_ensemblIDs.txt"))
adultback1 <- as.vector(adultback1$V1)
adultback2 <- read.table(paste0(anno, "werling2016_adultRep_ensemblIDs.txt"))
adultback2 <- as.vector(adultback2$V1)
prenatalback <- read.table(paste0(anno, "werling2016_prenatal_ensemblIDs.txt"))
prenatalback <- as.vector(prenatalback$V1)
# Adult1
adult1Female <- as.vector(read.table(paste0(anno, "werling2016_adult1_fem_p01.txt")))
adult1Female <- as.vector(adult1Female$V1)
adult1Male <- as.vector(read.table(paste0(anno, "werling2016_adult1_male_p01.txt")))
adult1Male <- as.vector(adult1Male$V1)

# Adult2
adult2Female <- read.table(paste0(anno, "werling2016_adultRep_fem_p01.txt"))
adult2Female <- as.vector(adult2Female$V1)
adult2Male <- read.table(paste0(anno, "werling2016_adultRep_male_p01.txt"))
adult2Male <- as.vector(adult2Male$V1)

# prenatal
prenatalFemale <- read.table(paste0(anno, "werling2016_prenatal_fem_p01.txt"))
prenatalFemale <- as.vector(prenatalFemale$V1)
prenatalMale <- read.table(paste0(anno,"werling2016_prenatal_male_p01.txt"))
prenatalMale<- as.vector(prenatalMale$V1)


# Read BrainSpan
anno <- "/Users/apple/Box/Li_Jiaxin/annotationLists/sexDex_brainSpan_windows-171025/"

# Background
BrainSpanback <- read.table(paste0(anno, "testedGenes_ensID_171025.txt"))
BrainSpanback <- as.vector(BrainSpanback$V1)

####
BrainSpanFemale <- read.table(paste0(anno, "femDex_ensID_Q0.05_171025.txt"))
BrainSpanFemale <- as.vector(BrainSpanFemale$V1)
BrainSpanMale <- read.table(paste0(anno, "maleDex_ensID_Q0.05_171025.txt"))
BrainSpanMale <- as.vector(BrainSpanMale$V1)

# Read Satterstrom 2020
anno <- "/Users/apple/Box/Li_Jiaxin/annotationLists/asdRisk_Satterstrom2020/"
SattBack <- read.table(paste0(anno, "ascGenesBackground-autosomal-ensIDs-2020.txt"))
SattBack <- as.vector(SattBack$V1)

Satt_ASDskew <- read.table(paste0(anno, "ascGenes-ASDskew-ensIDs-2020.txt"))
Satt_ASDskew <- as.vector(Satt_ASDskew$V1)
Satt_Cyto <- read.table(paste0(anno , "ascGenes-Cyto-ensIDs-2020.txt"))
Satt_Cyto <- as.vector(Satt_Cyto$V1)
Satt_DDIDskew <- read.table(paste0(anno, "ascGenes-DDIDskew-ensIDs-2020.txt"))
Satt_DDIDskew <- as.vector(Satt_DDIDskew$V1)
Satt_GER <- read.table(paste0(anno, "ascGenes-GER-ensIDs-2020.txt"))
Satt_GER <- as.vector(Satt_GER$V1)
Satt_NC <- read.table(paste0(anno, "ascGenes-NC-ensIDs-2020.txt"))
Satt_NC <- as.vector(Satt_NC$V1)
Satt_102 <- read.table(paste0(anno, "ascGenes102-ensIDs-2020.txt"))
Satt_102 <- as.vector(Satt_102$V1)

# Read Pari
anno <- "/mnt/sas0/AD/jli2274/annotationLists/asdModules_Parikshak2016/"

PariBack <- read.table(paste0(anno, "parikshak2016Background_ens.txt"))
PariBack <- as.vector(PariBack$V1)
Pari_m4 <- read.table(paste0(anno, "parikshak2016_ctx.m4.txt"))
Pari_m4 <- as.vector(Pari_m4$V1)
Pari_m9 <- read.table(paste0(anno, "parikshak2016_ctx.m9.txt"))
Pari_m9 <- as.vector(Pari_m9$V1)
Pari_m10 <- read.table(paste0(anno, "parikshak2016_ctx.m10.txt"))
Pari_m10 <- as.vector(Pari_m10$V1)
Pari_m16 <- read.table(paste0(anno, "parikshak2016_ctx.m16.txt"))
Pari_m16 <- as.vector(Pari_m16$V1)
Pari_m19 <- read.table(paste0(anno, "parikshak2016_ctx.m19.txt"))
Pari_m19 <- as.vector(Pari_m19$V1)
Pari_m20 <- read.table(paste0(anno, "parikshak2016_ctx.m20.txt"))
Pari_m20 <- as.vector(Pari_m20$V1)

# Read Velmeshev_2019
# Background
anno <- "/Users/apple/Box/Li_Jiaxin/annotationLists/Velmeshev_2019_DEGs.xls"
Velm <- read_excel(anno, sheet = 1 , na = "NA")
VelmCase <- as.vector(Velm[Velm$`Cell type` == celltype & Velm$`Fold change` < 0,]$`gene ID`)
VelmControl <- as.vector(Velm[Velm$`Cell type` == celltype & Velm$`Fold change` > 0,]$`gene ID`)
VelmBack <- as.vector(fdata$primerid)





