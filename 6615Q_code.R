### 1. Initial data manipulation ###

library("readxl")
library("ggplot2")
library("dplyr")
library("stringr")
library("devtools")
#devtools::install_github("uqrmaie1/admixtools")
library("admixtools")
library("genio")
library("data.table")
library("reshape2")
library("ggbreak")
library("RColorBrewer")
library("stringr")
library("tidyverse")
library("viridis")
library("ggrepel")
library("remotes")
library("extrafont")

# Load AADR meta files and divide the samples into batches of 500 years

AADR_meta <- read_excel("1240K_regions.xlsx")
colnames(AADR_meta)[7] <- "date.mean"
colnames(AADR_meta)[8] <- "date.std"
colnames(AADR_meta)[9] <- "full.date"


AADR_meta =  AADR_meta %>% 
  mutate(Period = case_when(date.mean >= 6500 ~ "6500+ BP",
                            date.mean >= 6000 & date.mean < 6500 ~ "6000 - 6499 BP",
                            date.mean >= 5500 & date.mean < 6000 ~ "5500 - 5999 BP",
                            date.mean >= 5000 & date.mean < 5500 ~ "5000 - 5499 BP",
                            date.mean >= 4500 & date.mean < 5000 ~ "4500 - 4999 BP",
                            date.mean >= 4000 & date.mean < 4500 ~ "4000 - 4499 BP",
                            date.mean >= 3500 & date.mean < 4000 ~ "3500 - 3999 BP",
                            date.mean >= 3000 & date.mean < 3500 ~ "3000 - 3499 BP",
                            date.mean >= 2500 & date.mean < 3000 ~ "2500 - 2999 BP",
                            date.mean >= 2000 & date.mean < 2500 ~ "2000 - 2499 BP",
                            date.mean >= 1500 & date.mean < 2000 ~ "1500 - 1999 BP",
                            date.mean >= 1000 & date.mean < 1500 ~ "1000 - 1499 BP",
                            date.mean >= 500 & date.mean < 1000 ~ "500 - 999 BP",
                            date.mean > 0 & date.mean < 500 ~ "1 - 499 BP",
                            date.mean == 0 ~ "Present"))

AADR_meta = AADR_meta[!(AADR_meta$Period=="6500+ BP"),]

AADR_meta$Period = factor(AADR_meta$Period, levels=c("Present", 
                                                   "1 - 499 BP", 
                                                   "500 - 999 BP", 
                                                   "1000 - 1499 BP",
                                                   "1500 - 1999 BP",
                                                   "2000 - 2499 BP",
                                                   "2500 - 2999 BP",
                                                   "3000 - 3499 BP",
                                                   "3500 - 3999 BP",
                                                   "4000 - 4499 BP",
                                                   "4500 - 4999 BP",
                                                   "5000 - 5499 BP",
                                                   "5500 - 5999 BP",
                                                   "6000 - 6499 BP"))



### 2.Choosing the populations/regions to be analysed ###

#Iberia
iber = subset(AADR_meta, Country=="Spain"|Country=="Portugal"|
                 Country=="Gibraltar")

iber_aDNA = iber[iber[,15] != "Present",]

mycolors_iber = colorRampPalette(brewer.pal(8, "BrBG"))(13)
mycolors_iber = rev(mycolors_iber)
                    
ggplot(iber_aDNA, aes(x=Period, fill=factor(..x..))) + 
  geom_histogram(stat="count", col="white") + 
  scale_fill_manual(values = mycolors_iber) + theme_bw() +
  theme(legend.position = "none") + 
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

#South East Asia
seasia = subset (AADR_meta, Country=="China"|Country=="Laos"|
                   Country=="Vietnam"|Country=="Thailand")

#Andes
andes = subset(AADR_meta, Country=="Peru"|Country=="Bolivia"|
                 Country=="Chile"& Lat. < (-22))

#Great Steppe
steppe = subset(AADR_meta, Country=="Russia"|Country=="Kazakhstan"|
                  Country=="Kyrgyzstan"|Country=="Mongolia"|
                  Country=="Ukraine")

steppe_aDNA = steppe[steppe[,15] != "Present",]

mycolors_steppe = colorRampPalette(brewer.pal(8, "YlGnBu"))(13)
mycolors_steppe = rev(mycolors_steppe)

ggplot(steppe_aDNA, aes(x=Period, fill=factor(..x..))) + 
  geom_histogram(stat="count", col="white") + 
  scale_fill_viridis(discrete=TRUE) + theme_bw() +
  theme(legend.position = "none") + 
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

#East Africa
eafrica = subset(AADR_meta, Country=="Malawi"|Country=="Tanzania"|
                  Country=="Kenya"|Country=="Ethiopia")

#NOTE: The commands below apply only to the Iberian cline. Nevertheless, 
# workflow was exactly the same for all the clines listed above. 


# Create separate subsets of 500 years for the specific populations
iber0 = subset(iber, Period =="Present")
iber500 = subset(iber, Period =="1 - 499 BP")
iber1000 = subset(iber, Period =="500 - 999 BP")
iber1500 = subset(iber, Period =="1000 - 1499 BP")
iber2000 = subset(iber, Period =="1500 - 1999 BP")
iber2500 = subset(iber, Period =="2000 - 2499 BP")
iber3000 = subset(iber, Period =="2500 - 2999 BP")
iber3500 = subset(iber, Period =="3000 - 3499 BP")
iber4000 = subset(iber, Period =="3500 - 3999 BP")
iber4500 = subset(iber, Period =="4000 - 4499 BP")
iber5000 = subset(iber, Period =="4500 - 4999 BP")
iber5500 = subset(iber, Period =="5000 - 5499 BP")
iber6000 = subset(iber, Period =="5500 - 5999 BP")
iber6500 = subset(iber, Period =="6000 - 6499 BP")

names(iber0)<-str_replace_all(names(iber0), c(" " = "." , "," = "" ))
names(iber500)<-str_replace_all(names(iber500), c(" " = "." , "," = "" ))
names(iber1000)<-str_replace_all(names(iber1000), c(" " = "." , "," = "" ))
names(iber1500)<-str_replace_all(names(iber1500), c(" " = "." , "," = "" ))
names(iber2000)<-str_replace_all(names(iber2000), c(" " = "." , "," = "" ))
names(iber2500)<-str_replace_all(names(iber2500), c(" " = "." , "," = "" ))
names(iber3000)<-str_replace_all(names(iber3000), c(" " = "." , "," = "" ))
names(iber3500)<-str_replace_all(names(iber3500), c(" " = "." , "," = "" ))
names(iber4000)<-str_replace_all(names(iber4000), c(" " = "." , "," = "" ))
names(iber4500)<-str_replace_all(names(iber4500), c(" " = "." , "," = "" ))
names(iber5000)<-str_replace_all(names(iber5000), c(" " = "." , "," = "" ))
names(iber5500)<-str_replace_all(names(iber5500), c(" " = "." , "," = "" ))
names(iber6000)<-str_replace_all(names(iber6000), c(" " = "." , "," = "" ))
names(iber6500)<-str_replace_all(names(iber6500), c(" " = "." , "," = "" ))


# Export the lists of individuals
iber0_ind = unique(iber0$Version.ID)
iber500_ind = unique(iber500$Version.ID)
iber1000_ind = unique(iber1000$Version.ID)
iber1500_ind = unique(iber1500$Version.ID)
iber2000_ind = unique(iber2000$Version.ID)
iber2500_ind = unique(iber2500$Version.ID)
iber3000_ind = unique(iber3000$Version.ID)
iber3500_ind = unique(iber3500$Version.ID)
iber4000_ind = unique(iber4000$Version.ID)
iber4500_ind = unique(iber4500$Version.ID)
iber5000_ind = unique(iber5000$Version.ID)
iber5500_ind = unique(iber5500$Version.ID)
iber6000_ind = unique(iber6000$Version.ID)
iber6500_ind = unique(iber6500$Version.ID)


### 3. Conversion to PLINK format ###

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline0",
  inds = iber0_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline500",
  inds = iber500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline1000",
  inds = iber1000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline1500",
  inds = iber1500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline2000",
  inds = iber2000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline2500",
  inds = iber2500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline3000",
  inds = iber3000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline3500",
  inds = iber3500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline4000",
  inds = iber4000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline4500",
  inds = iber4500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline5000",
  inds = iber5000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline5500",
  inds = iber5500_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline6000",
  inds = iber6000_ind,
  pops = NULL,
  verbose = TRUE)

packedancestrymap_to_plink(
  inpref = "v44.3_1240K_public",
  outpref = "iberian_cline6500",
  inds = iber6500_ind,
  pops = NULL,
  verbose = TRUE)

# See 6615Q_PLINK_conversion.txt for PLINK workflow.


### 4. Generating inputs for Schreiber's algorithm ###

# 4.1 prepare frequency files

ib0_frq = read.table("iberian_cline0.frq", header = TRUE)
ib0_frq = ib0_frq[ , c(1,2,5,6)]
ib0_frq = ib0_frq[complete.cases(ib0_frq), ]
Count = ib0_frq$MAF * ib0_frq$NCHROBS
Count = round(Count, digits = 0)
ib0_frq = cbind(ib0_frq, Count)
ib0_frq = ib0_frq[ , c(1,2,5,4)]

ib500_frq = read.table("iberian_cline500.frq", header = TRUE)
ib500_frq = ib500_frq[ , c(1,2,5,6)]
ib500_frq = ib500_frq[complete.cases(ib500_frq), ]
Count = ib500_frq$MAF * ib500_frq$NCHROBS
Count = round(Count, digits = 0)
ib500_frq = cbind(ib500_frq, Count)
ib500_frq = ib500_frq[ , c(1,2,5,4)]

ib1000_frq = read.table("iberian_cline1000.frq", header = TRUE)
ib1000_frq = ib1000_frq[ , c(1,2,5,6)]
ib1000_frq = ib1000_frq[complete.cases(ib1000_frq), ]
Count = ib1000_frq$MAF * ib1000_frq$NCHROBS
Count = round(Count, digits = 0)
ib1000_frq = cbind(ib1000_frq, Count)
ib1000_frq = ib1000_frq[ , c(1,2,5,4)]

ib1500_frq = read.table("iberian_cline1500.frq", header = TRUE)
ib1500_frq = ib1500_frq[ , c(1,2,5,6)]
ib1500_frq = ib1500_frq[complete.cases(ib1500_frq), ]
Count = ib1500_frq$MAF * ib1500_frq$NCHROBS
Count = round(Count, digits = 0)
ib1500_frq = cbind(ib1500_frq, Count)
ib1500_frq = ib1500_frq[ , c(1,2,5,4)]

ib2000_frq = read.table("iberian_cline2000.frq", header = TRUE)
ib2000_frq = ib2000_frq[ , c(1,2,5,6)]
ib2000_frq = ib2000_frq[complete.cases(ib2000_frq), ]
Count = ib2000_frq$MAF * ib2000_frq$NCHROBS
Count = round(Count, digits = 0)
ib2000_frq = cbind(ib2000_frq, Count)
ib2000_frq = ib2000_frq[ , c(1,2,5,4)]

ib2500_frq = read.table("iberian_cline2500.frq", header = TRUE)
ib2500_frq = ib2500_frq[ , c(1,2,5,6)]
ib2500_frq = ib2500_frq[complete.cases(ib2500_frq), ]
Count = ib2500_frq$MAF * ib2500_frq$NCHROBS
Count = round(Count, digits = 0)
ib2500_frq = cbind(ib2500_frq, Count)
ib2500_frq = ib2500_frq[ , c(1,2,5,4)]

ib3000_frq = read.table("iberian_cline3000.frq", header = TRUE)
ib3000_frq = ib3000_frq[ , c(1,2,5,6)]
ib3000_frq = ib3000_frq[complete.cases(ib3000_frq), ]
Count = ib3000_frq$MAF * ib3000_frq$NCHROBS
Count = round(Count, digits = 0)
ib3000_frq = cbind(ib3000_frq, Count)
ib3000_frq = ib3000_frq[ , c(1,2,5,4)]

ib3500_frq = read.table("iberian_cline3500.frq", header = TRUE)
ib3500_frq = ib3500_frq[ , c(1,2,5,6)]
ib3500_frq = ib3500_frq[complete.cases(ib3500_frq), ]
Count = ib3500_frq$MAF * ib3500_frq$NCHROBS
Count = round(Count, digits = 0)
ib3500_frq = cbind(ib3500_frq, Count)
ib3500_frq = ib3500_frq[ , c(1,2,5,4)]

ib4000_frq = read.table("iberian_cline4000.frq", header = TRUE)
ib4000_frq = ib4000_frq[ , c(1,2,5,6)]
ib4000_frq = ib4000_frq[complete.cases(ib4000_frq), ]
Count = ib4000_frq$MAF * ib4000_frq$NCHROBS
Count = round(Count, digits = 0)
ib4000_frq = cbind(ib4000_frq, Count)
ib4000_frq = ib4000_frq[ , c(1,2,5,4)]

ib4500_frq = read.table("iberian_cline4500.frq", header = TRUE)
ib4500_frq = ib4500_frq[ , c(1,2,5,6)]
ib4500_frq = ib4500_frq[complete.cases(ib4500_frq), ]
Count = ib4500_frq$MAF * ib4500_frq$NCHROBS
Count = round(Count, digits = 0)
ib4500_frq = cbind(ib4500_frq, Count)
ib4500_frq = ib4500_frq[ , c(1,2,5,4)]

ib5000_frq = read.table("iberian_cline5000.frq", header = TRUE)
ib5000_frq = ib5000_frq[ , c(1,2,5,6)]
ib5000_frq = ib5000_frq[complete.cases(ib5000_frq), ]
Count = ib5000_frq$MAF * ib5000_frq$NCHROBS
Count = round(Count, digits = 0)
ib5000_frq = cbind(ib5000_frq, Count)
ib5000_frq = ib5000_frq[ , c(1,2,5,4)]

ib5500_frq = read.table("iberian_cline5500.frq", header = TRUE)
ib5500_frq = ib5500_frq[ , c(1,2,5,6)]
ib5500_frq = ib5500_frq[complete.cases(ib5500_frq), ]
Count = ib5500_frq$MAF * ib5500_frq$NCHROBS
Count = round(Count, digits = 0)
ib5500_frq = cbind(ib5500_frq, Count)
ib5500_frq = ib5500_frq[ , c(1,2,5,4)]

ib6000_frq = read.table("iberian_cline6000.frq", header = TRUE)
ib6000_frq = ib6000_frq[ , c(1,2,5,6)]
ib6000_frq = ib6000_frq[complete.cases(ib6000_frq), ]
Count = ib6000_frq$MAF * ib6000_frq$NCHROBS
Count = round(Count, digits = 0)
ib6000_frq = cbind(ib6000_frq, Count)
ib6000_frq = ib6000_frq[ , c(1,2,5,4)]

# 4.2 Standardize the sample - assure that every set contains the same SNPs.

ib0_frq = ib0_frq[ib0_frq$SNP %in% ib500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib1000_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib1500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib2000_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib2500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib3000_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib3500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib4000_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib4500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib5000_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib5500_frq$SNP,]
ib0_frq = ib0_frq[ib0_frq$SNP %in% ib6000_frq$SNP,]

ib500_frq = ib500_frq[ib500_frq$SNP %in% ib0_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib1000_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib1500_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib2000_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib2500_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib3000_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib3500_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib4000_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib4500_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib5000_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib5500_frq$SNP,]
ib500_frq = ib500_frq[ib500_frq$SNP %in% ib6000_frq$SNP,]

ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib0_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib1500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib2000_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib2500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib3000_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib3500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib4000_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib4500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib5000_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib5500_frq$SNP,]
ib1000_frq = ib1000_frq[ib1000_frq$SNP %in% ib6000_frq$SNP,]

ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib0_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib1000_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib500_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib2000_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib2500_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib3000_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib3500_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib4000_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib4500_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib5000_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib5500_frq$SNP,]
ib1500_frq = ib1500_frq[ib1500_frq$SNP %in% ib6000_frq$SNP,]

ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib0_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib1500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib1000_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib2500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib3000_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib3500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib4000_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib4500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib5000_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib5500_frq$SNP,]
ib2000_frq = ib2000_frq[ib2000_frq$SNP %in% ib6000_frq$SNP,]

ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib0_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib1000_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib500_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib2000_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib1500_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib3000_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib3500_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib4000_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib4500_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib5000_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib5500_frq$SNP,]
ib2500_frq = ib2500_frq[ib2500_frq$SNP %in% ib6000_frq$SNP,]

ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib0_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib1500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib1000_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib2500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib2000_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib3500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib4000_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib4500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib5000_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib5500_frq$SNP,]
ib3000_frq = ib3000_frq[ib3000_frq$SNP %in% ib6000_frq$SNP,]

ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib0_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib1000_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib500_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib2000_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib1500_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib3000_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib2500_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib4000_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib4500_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib5000_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib5500_frq$SNP,]
ib3500_frq = ib3500_frq[ib3500_frq$SNP %in% ib6000_frq$SNP,]

ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib0_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib1500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib1000_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib2500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib2000_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib3500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib3000_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib4500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib5000_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib5500_frq$SNP,]
ib4000_frq = ib4000_frq[ib4000_frq$SNP %in% ib6000_frq$SNP,]

ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib0_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib1000_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib500_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib2000_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib1500_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib3000_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib2500_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib4000_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib3500_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib5000_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib5500_frq$SNP,]
ib4500_frq = ib4500_frq[ib4500_frq$SNP %in% ib6000_frq$SNP,]

ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib0_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib1500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib1000_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib2500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib2000_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib3500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib3000_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib4500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib4000_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib5500_frq$SNP,]
ib5000_frq = ib5000_frq[ib5000_frq$SNP %in% ib6000_frq$SNP,]

ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib0_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib1000_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib500_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib2000_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib1500_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib3000_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib2500_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib4000_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib3500_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib5000_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib4500_frq$SNP,]
ib5500_frq = ib5500_frq[ib5500_frq$SNP %in% ib6000_frq$SNP,]

ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib0_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib1500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib1000_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib2500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib2000_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib3500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib3000_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib4500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib4000_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib5500_frq$SNP,]
ib6000_frq = ib6000_frq[ib6000_frq$SNP %in% ib5000_frq$SNP,]


# 4.3 Write down the start and end points for each time period (based on meta 
# data provided by AADR for the analysed population).

starting_times = c(-5500, -4906, -4240, -3935, -3461, -3013, 
                   -2471, -2045, -1550,  -971,  -492,  -400, 0)

ending_times = c(-6123, -5535, -5051, -4507, -4001, -3517, 
                 -3044, -2527, -1947, -1504, -1094,  -400, 0)


# 4.4 Create a function to generate the outputs. Input the number of row from
# the data sets generated above.

iber_files_input = function(i) {
  result = rbind(ib6000_frq[i,], ib5500_frq[i,], ib5000_frq[i,], ib4500_frq[i,],
                 ib4000_frq[i,], ib3500_frq[i,], ib3000_frq[i,], ib2500_frq[i,], 
                 ib2000_frq[i,], ib1500_frq[i,], ib1000_frq[i,], ib500_frq[i,],
                 ib0_frq[i,])
  df = data.frame(result)
  df = df[, c(3,4)]
  df = cbind(df, ending_times, starting_times)
  write.table(df, paste0(ib500_frq[i,2], ".txt"), col.names = FALSE, row.names = FALSE, sep="\t")
  return(df)
}


# 4.5 Loop the sr_input function over desired number of rows. 

for (i in 1:40) {
  iber_files_input(i)
}


# 4.5a Generating one file for the enitre cline/ sample


iber_one_file = function(i) {
  a_freq = cbind(ib6000_frq[i,], ib5500_frq[i,c(3,4)], ib5000_frq[i,c(3,4)], ib4500_frq[i,c(3,4)],
                 ib4000_frq[i,c(3,4)], ib3500_frq[i,c(3,4)], ib3000_frq[i,c(3,4)], ib2500_frq[i,c(3,4)], 
                 ib2000_frq[i,c(3,4)], ib1500_frq[i,c(3,4)], ib1000_frq[i,c(3,4)], ib500_frq[i,c(3,4)],
                 ib0_frq[i,c(3,4)])
  df = data.frame(a_freq)
  write.table(df)
  return(df)
}


iber_input = sapply(1:881482, iber_one_file)
iber_input = t(iber_input)
iber_input = as.data.frame(iber_input)
colnames(iber_input) <- c("CHR", "SNP", "Count6000", "N6000","Count5500", "N5500",
                     "Count5000", "N5000","Count4500", "N4500","Count4000", "N4000",
                     "Count3500", "N3500","Count3000", "N3000","Count2500", "N2500",
                     "Count2000", "N2000","Count1500", "N1500","Count1000", "N1000",
                     "Count500", "N500", "Count0", "N0")

iber_input = as.data.table(iber_input)
fwrite(iber_input, "iber_non_pruned.csv")
iber_non_pruned = read.csv(file = "iber_non_pruned.csv")

iber_pruned = read.delim("iber_pruned.prune.in", header = FALSE)
iber_chr6 = iber_chr6[iber_chr6$SNP %in% iber_pruned$V1,]

# Pruning (from output file generated with PLINK LD pruning function)

iber_pruned = read.delim("iber_pruned.prune.in", header = FALSE)
iber_input = iber_non_pruned[iber_non_pruned$SNP %in% iber_pruned$V1,]

fwrite(iber_input, "iber_pruned_input.csv")
iber_input = read.csv(file = "iber_input.csv")





### 5. Generating frequency bins for the analysis ###

iber_modern_freq = read.table("iberian_cline0.frq", header = TRUE)
iber_modern_freq = iber_modern_freq[ , c(1,2,5)]
iber_modern_freq = iber_modern_freq[complete.cases(iber_modern_freq), ]
iber_modern_freq = iber_modern_freq[iber_modern_freq$SNP %in% iber_count$SNP,]

iber_modern_freq =  iber_modern_freq %>% 
  mutate(freq_bin= case_when(MAF >= 0.95 & 1 >= MAF ~ "95-100%",
                             MAF >= 0.9 & MAF < 0.95 ~ "90-95%",
                             MAF >= 0.85 & MAF < 0.9 ~ "85-90%",
                             MAF >= 0.8 & MAF < 0.85 ~ "80-85%",
                             MAF >= 0.75 & MAF < 0.8 ~ "75-80%",
                             MAF >= 0.7 & MAF < 0.75 ~ "70-75%",
                             MAF >= 0.65 & MAF < 0.7 ~ "65-70%",
                             MAF >= 0.6 & MAF < 0.65 ~ "60-65%",
                             MAF >= 0.55 & MAF < 0.6 ~ "55-60%",
                             MAF >= 0.5 & MAF < 0.55 ~ "50-55%",
                             MAF >= 0.45 & MAF < 0.5 ~ "45-50%",
                             MAF >= 0.4 & MAF < 0.45 ~ "40-45%",
                             MAF >= 0.35 & MAF < 0.4 ~ "35-40%",
                             MAF >= 0.3 & MAF < 0.35 ~ "30-35%",
                             MAF >= 0.25 & MAF < 0.3 ~ "25-30%",
                             MAF >= 0.2 & MAF < 0.25 ~ "20-25%",
                             MAF >= 0.15 & MAF < 0.2 ~ "15-20%",
                             MAF >= 0.1 & MAF < 0.15 ~ "10-15%",
                             MAF >= 0.05 & MAF< 0.1 ~ "5-10%",
                             MAF >= 0 & MAF < 0.05 ~ "0-5%"))

ib_bin1 = subset(iber_modern_freq, freq_bin=="0-5%")
ib_bin2 = subset(iber_modern_freq, freq_bin=="5-10%")
ib_bin3 = subset(iber_modern_freq, freq_bin=="10-15%")
ib_bin4 = subset(iber_modern_freq, freq_bin=="15-20%")
ib_bin5 = subset(iber_modern_freq, freq_bin=="20-25%")
ib_bin6 = subset(iber_modern_freq, freq_bin=="25-30%")
ib_bin7 = subset(iber_modern_freq, freq_bin=="30-35%")
ib_bin8 = subset(iber_modern_freq, freq_bin=="35-40%")
ib_bin9 = subset(iber_modern_freq, freq_bin=="40-45%")
ib_bin10 = subset(iber_modern_freq, freq_bin=="45-50%")
ib_bin11 = subset(iber_modern_freq, freq_bin=="50-55%")
ib_bin12 = subset(iber_modern_freq, freq_bin=="55-60%")
ib_bin13 = subset(iber_modern_freq, freq_bin=="60-65%")
ib_bin14 = subset(iber_modern_freq, freq_bin=="65-70%")
ib_bin15 = subset(iber_modern_freq, freq_bin=="70-75%")
ib_bin16 = subset(iber_modern_freq, freq_bin=="75-80%")
ib_bin17 = subset(iber_modern_freq, freq_bin=="80-85%")
ib_bin18 = subset(iber_modern_freq, freq_bin=="85-90%")
ib_bin19 = subset(iber_modern_freq, freq_bin=="90-95%")
ib_bin20 = subset(iber_modern_freq, freq_bin=="95-100%")



# Sample 100 random rows from each bin

ib_bin_list = list(ib_bin1, ib_bin2, ib_bin3, ib_bin4, ib_bin5, ib_bin6, ib_bin7, 
                ib_bin8, ib_bin9, ib_bin10, ib_bin11, ib_bin12, ib_bin13, 
                ib_bin14, ib_bin15, ib_bin16, ib_bin17, ib_bin18, 
                ib_bin19, ib_bin20)

sample_bins = function(i) {
  obj_name = deparse(substitute(i))
  i = i[sample(nrow(i), 100), ]
  i = iber_input[iber_input$SNP %in% i$SNP,]
  write.csv(i, file = paste0(obj_name),  row.names = FALSE)
}



### Create sets of SNPs

snp_position = read.table("v44.3_1240K_public.snp")

HLA_pos = subset(snp_position, V4 >= 29677984 & 33485635 >= V4)
iber_HLA = iber_count[iber_count$SNP %in% HLA_pos$V1,]
fwrite(iber_HLA, "iber_HLA")

sel_cand_pos = subset(snp_position, V4 >= 169367970 & 169396575 >= V4)

fwrite(AADR_meta, "AADR_locations.csv")


### Analyse trajectories for candidate selection SNPs ###

sel_SNPs = c("rs4988235", "rs1426654", "rs3827760", "rs17822931", "rs1871534",
             "rs1050828", "rs1050829", "rs5030868", "rs137852314", "rs76723693",
             "rs137852328", "rs372091", "rs33930165", "rs33950507", "rs334",
             "rs2814778", "rs1050501", "rs8177374", "rs4951074", "rs10900585",
             "rs2230345", "rs1800890", "rs2334880", "rs8176719", "rs8176746",
             "rs3092945", "rs201346212")


candidate_SNPs = subset(count, SNP %in% sel_SNPs)


### See the respective files for each cline for output analysis ###


