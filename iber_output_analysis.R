
### 5. Generating frequency bins for the analysis ###

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/Inputs/pruned_inputs")
iber_input = read.table("iber_pruned_input.csv", header = TRUE, sep=",")

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
iber_modern_freq = read.table("iberian_cline0.frq", header = TRUE)

iber_modern_freq = iber_modern_freq[ , c(1,2,5)]
iber_modern_freq = iber_modern_freq[complete.cases(iber_modern_freq), ]
iber_modern_freq = iber_modern_freq[iber_modern_freq$SNP %in% iber_input$SNP,]

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

# Sample 100 random rows from each bin
iber_neutrality_snps = iber_modern_freq %>% group_by(freq_bin) %>% slice_sample(n=100)
iber_neutrality_input = iber_input[iber_input$SNP %in% iber_neutrality_snps$SNP, ]

setwd("C:/Users/ACER/Desktop")
fwrite(iber_neutrality_input, "iber_neutrality_input.csv")


#Read the results for bins' analysis
iber_neutrality_results = read.table("iber_neutrality_results.csv", header = FALSE, sep=",")
colnames(iber_neutrality_results) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")

iber_neutrality_snps = iber_modern_freq[,c(1,2,4)]
iber_neutrality_results = merge(iber_neutrality_results, iber_modern_freq, by="SNP")

#Check and plot the distributions of alpha1 and alpha2 
setDT(iber_neutrality_results)
iber_neutr_a1 = as.data.frame(iber_neutrality_results[, as.list(summary(alpha1_mean)), by = freq_bin])
iber_neutr_a2 = as.data.frame(iber_neutrality_results[, as.list(summary(alpha2_mean)), by = freq_bin])


iber_neutr_a1$freq_bin = factor(iber_neutr_a1$freq_bin,
                                levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                         "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                         "85-90%","90-95%","95-100%"))

iber_neutr_a2$freq_bin = factor(iber_neutr_a2$freq_bin,
                                levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                         "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                         "85-90%","90-95%","95-100%"))



#alpha1 mean in different bins barplot
ggplot(data=iber_neutr_a1, aes(x=freq_bin, y=Median)) + geom_bar(stat="identity")

#alpha2 mean in different bins barplot
ggplot(data=iber_neutr_a2, aes(x=freq_bin, y=Median)) + geom_bar(stat="identity")

#Plot distribution of alpha1 and alpha2 within the bins
iber_neutr_a1$Coefficient = "Alpha 1" 
iber_neutr_a2$Coefficient = "Alpha 2"
iber_neutr_merged = rbindlist(list(iber_neutr_a1,iber_neutr_a2))

ggplot(data = iber_neutr_merged, aes(x = freq_bin, y = Median, color = Coefficient, fill = Coefficient)) + 
  geom_bar(stat = "identity", alpha=1, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = c("#f6b60d","#00AFBB")) + 
  xlab("Frequency bin") + theme_bw() +
  scale_x_discrete(breaks = levels(iber_neutr_merged$freq_bin)[c(T, rep(F, 8))]) +
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"), 
        panel.border = element_blank()) +
labs(title="Iberian cline: median of selection coefficients in different 'neutrality bins'") +
  labs(subtitle="n = 2000 SNPs from 'neutrality bins' (20 bins of 100 SNP sampled by frequency every 5%)") +
  scale_fill_manual(name="Coefficient",labels=c(expression(" "*alpha[1]), expression("  "*alpha[2])), values=c("#E7B800","#00AFBB"))


#Plot alpha1 vs alpha2
setwd("C:/Users/ACER/Desktop")
iber_neutrality_results = read.table("iber_neutrality_results.csv", header = FALSE, sep=",")
colnames(iber_neutrality_results) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")
iber_neutrality_results = iber_neutrality_results[-c(1403,137,133,851,827,311,31, 1057,415), ]

iber_neutrality_results$Selection[iber_neutrality_results$alpha2_mean > iber_neutrality_results$alpha1_mean & iber_neutrality_results$alpha1_mean > 0] = "Directional"
iber_neutrality_results$Selection[iber_neutrality_results$alpha1_mean > iber_neutrality_results$alpha2_mean & iber_neutrality_results$alpha1_mean > 0] = "Overdominant"
iber_neutrality_results$Selection[iber_neutrality_results$alpha2_mean > iber_neutrality_results$alpha1_mean & 0 > iber_neutrality_results$alpha1_mean] = "Underdominant"

table(iber_neutrality_results$Selection)

ggplot(data=iber_neutrality_results, aes(x=alpha1_mean, y=alpha2_mean, color=Selection)) + 
  geom_point(alpha=0.5) + geom_smooth(method=lm , color="#404040", se=FALSE, linetype="dashed") + 
  theme_bw() + xlab(expression(alpha[1]*" coefficient")) + ylab(expression(alpha[2]*" coefficient")) +
  scale_color_manual(values=c("#00AFBB", "#f6b60d")) +
  labs(title="Iberian cline: correlation between selection coefficients") +
  labs(subtitle="n = 2000 SNPs from 'neutrality bins' (20 bins of 100 SNP sampled by frequency every 5%)") +
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"), 
        panel.border = element_blank()) + ylim(-50,1600)



#Read the results for the entire chromosome 6
setwd("C:/Users/ACER/Desktop")
iber_chr6_res = read.table("iber_chr6_results.csv", header = FALSE, sep=",")
colnames(iber_chr6_res) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")



setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
snp_position = read.table("iberian_cline0.map")
snp_position = snp_position[, c(1:4)]
colnames(snp_position) = c("CHR","SNP","gpos", "phpos")
iber_chr6_snp = snp_position[snp_position$SNP %in% iber_chr6_res$SNP, ]
iber_chr6_res = merge(iber_chr6_snp, iber_chr6_res, by="SNP")

iber_modern_freq = iber_modern_freq[,c(2,4)]
iber_chr6_res = merge(iber_chr6_res, iber_modern_freq, by="SNP")

iber_chr6_res$freq_bin = factor(iber_chr6_res$freq_bin,
                                levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                         "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                         "85-90%","90-95%","95-100%"))

#Look for outliers in the entire data set

#Calculate IQR
iber_neutr_a1$IQR = iber_neutr_a1$`3rd Qu.` - iber_neutr_a1$`1st Qu.`
iber_neutr_a2$IQR = iber_neutr_a2$`3rd Qu.` - iber_neutr_a2$`1st Qu.`

#Calculate upper fence for strong outliers
iber_neutr_a1$UFSO = iber_neutr_a1$`3rd Qu.` + 3 * iber_neutr_a1$IQR
iber_neutr_a2$UFSO = iber_neutr_a2$`3rd Qu.` + 3 * iber_neutr_a2$IQR

#Calculate lower fence for strong outliers
iber_neutr_a1$LFSO = iber_neutr_a1$`1st Qu.` - 3 * iber_neutr_a1$IQR
iber_neutr_a2$LFSO = iber_neutr_a2$`1st Qu.` - 3 * iber_neutr_a2$IQR

#UF for a1
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="0-5%"] = 46.59410
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="5-10%"] = 43.76438
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="10-15%"] = 75.57841 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="15-20%"] = 64.40388 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="20-25%"] = 70.00471 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="25-30%"] = 66.30027 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="30-35%"] = 83.23360
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="35-40%"] = 108.45748
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="40-45%"] = 111.47388
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="45-50%"] = 109.17985 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="50-55%"] = 136.85276 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="55-60%"] = 100.74769 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="60-65%"] = 1061.77878 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="65-70%"] = 960.45690 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="70-75%"] = 1488.68735
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="75-80%"] = 1975.63059 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="80-85%"] = 1737.92052 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="85-90%"] = 1746.76125 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="90-95%"] = 1668.60258 
iber_chr6_res$UFSO.a1[iber_chr6_res$freq_bin=="95-100%"] = 1584.01352 

#UF for a2
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="0-5%"] = 42.37011
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="5-10%"] = 42.76671
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="10-15%"] = 49.75827
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="15-20%"] = 36.15449
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="20-25%"] = 45.51204
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="25-30%"] = 34.11629
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="30-35%"] = 48.43919
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="35-40%"] = 45.80850
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="40-45%"] = 52.47080
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="45-50%"] = 42.12496
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="50-55%"] = 93.28530
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="55-60%"] = 88.31204
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="60-65%"] = 985.86676
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="65-70%"] = 925.73076
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="70-75%"] = 1374.06552
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="75-80%"] = 1822.77445
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="80-85%"] = 1595.78546
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="85-90%"] = 1620.75161
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="90-95%"] = 1577.20087
iber_chr6_res$UFSO.a2[iber_chr6_res$freq_bin=="95-100%"] = 1490.73485

#Find & classify a1 outliers

iber_chr6_res$outa1[iber_chr6_res$UFSO.a1 > iber_chr6_res$alpha1_mean | iber_chr6_res$alpha1_mean < iber_chr6_res$alpha2_mean] = "Regular SNP"
iber_chr6_res$outa1[iber_chr6_res$UFSO.a1 < iber_chr6_res$alpha1_mean & iber_chr6_res$alpha1_mean > iber_chr6_res$alpha2_mean] = "Candidate SNP"

table(iber_chr6_res$outa1)

d = subset(iber_chr6_res, outa1=="Candidate SNP" & alpha1_mean > alpha2_mean)

#Find & classify a2 outliers

iber_chr6_res$outa2[iber_chr6_res$UFSO.a2 > iber_chr6_res$alpha2_mean | iber_chr6_res$alpha1_mean > iber_chr6_res$alpha2_mean] = "Regular SNP"
iber_chr6_res$outa2[iber_chr6_res$UFSO.a2 < iber_chr6_res$alpha2_mean & iber_chr6_res$alpha1_mean < iber_chr6_res$alpha2_mean] = "Candidate SNP"

table(iber_chr6_res$outa2)

e = subset(iber_chr6_res, outa2=="Candidate SNP" & alpha1_mean < alpha2_mean)

#Estimate the amount of selection candidates within and outside HLA

hla_out_a1 = subset(d, phpos >= 29677984 & 33485635 >= phpos)
hla_non_out_a1 = subset(d, phpos < 29677984 | 33485635 < phpos)

hla_out_a2 = subset(e, phpos >= 29677984 & 33485635 >= phpos)
hla_non_out_a2 = subset(e, phpos < 29677984 | 33485635 < phpos)

hla_SNPs = subset(iber_chr6_res, phpos >= 29677984 & 33485635 >= phpos)
non_hla_SNPs = subset(iber_chr6_res, phpos < 29677984 | 33485635 < phpos)

#Classify selection type

iber_chr6_res$sel[iber_chr6_res$alpha2_mean > iber_chr6_res$alpha1_mean & iber_chr6_res$alpha1_mean > 0] = "Direct"
iber_chr6_res$sel[iber_chr6_res$alpha1_mean > iber_chr6_res$alpha2_mean & iber_chr6_res$alpha1_mean > 0] = "Overdominant"
iber_chr6_res$sel[iber_chr6_res$alpha2_mean > iber_chr6_res$alpha1_mean & 0 > iber_chr6_res$alpha1_mean] = "Underdominant"

table(iber_chr6_res$sel)


# Find top hits for alpha1
iber_candidates_a1 = subset(iber_chr6_res, outa1=="Candidate SNP"& alpha1_mean > alpha2_mean)
iber_candidates_a1$diff_a1 = iber_candidates_a1$alpha1_mean/iber_candidates_a1$UFSO.a1
iber_candidates_a1 = iber_candidates_a1[order(iber_candidates_a1$diff_a1, decreasing = TRUE), ] 

iber_top10_a1 = top_n(iber_candidates_a1, 10, diff_a1)
iber_chr6_res$outa1[iber_chr6_res$SNP %in% iber_top10_a1$SNP] = "Top 10 hits"

  #Export the results
#iber_top10_export = iber_top10_a1[,c(1,4,5,14)]
#fwrite(iber_top10_export, "iber_top10.csv")

hla_top_iber_a1 = subset(iber_top10_a1, phpos >= 29677984 & 33485635 >= phpos)

ggplot(iber_chr6_res, aes(x=phpos, y=alpha1_mean, color=outa1)) + 
  geom_point(alpha=0.15) +
  geom_vline(xintercept = 29677984, alpha=0.2, linetype="longdash") +
  geom_vline(xintercept = 33485635, alpha=0.2, linetype="longdash") +
  geom_point(data = subset(iber_chr6_res, outa1 == "Candidate SNP"), alpha=0.6) + 
  geom_point(data = subset(iber_chr6_res, outa1 == "Top 10 hits"), alpha=0.8) +
  geom_rug(data = subset(iber_chr6_res, outa1 == "Candidate SNP"), sides="b") +
  geom_rug(data = subset(iber_chr6_res, outa1 == "Top 10 hits"), sides="b") +
  labs(title="Iberian cline: candidate SNPs for strong selection on heterozygotes",
       subtitle="Chr. 6 selection scan conducted on 10 365 SNPs in 109 modern DNA and 393 aDNA samples") +
  xlab("Physical position (in BP units)") + ylab(expression(alpha[1]*' coefficient')) +
  scale_color_manual(values=c("#47abd8","#D3D3D3","#D01B1B")) + 
  theme_bw() + theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) +
  geom_label_repel(data=iber_chr6_res[iber_chr6_res$outa1 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, 
                   xlim  = c(0, 19677984, 33485635, Inf), family = "Bahnschrift", show.legend=FALSE) + 
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  labs(color="SNP type") + 
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"))


table(iber_candidates_a1$sel)


# Find top hits for alpha2
iber_candidates_a2 = subset(iber_chr6_res, outa2=="Candidate SNP" & alpha2_mean > alpha1_mean)
iber_candidates_a2$diff = iber_candidates_a2$alpha2_mean/iber_candidates_a2$UFSO.a2
iber_candidates_a2 = iber_candidates_a2[order(iber_candidates_a2$diff, decreasing = TRUE), ] 

iber_top10_a2 = top_n(iber_candidates_a2, 10, diff)
iber_chr6_res$outa2[iber_chr6_res$SNP %in% iber_top10_a2$SNP] = "Top 10 hits"

hla_top_iber_a2 = subset(iber_top10_a2, phpos >= 29677984 & 33485635 >= phpos)

ggplot(iber_chr6_res, aes(x=phpos, y=alpha2_mean, color=outa2)) + 
  geom_point(alpha=0.15) +
  geom_vline(xintercept = 29677984, alpha=0.2, linetype="longdash") +
  geom_vline(xintercept = 33485635, alpha=0.2, linetype="longdash") +
  geom_point(data = subset(iber_chr6_res, outa2 == "Candidate SNP"), alpha=0.6) + 
  geom_point(data = subset(iber_chr6_res, outa2 == "Top 10 hits"), alpha=0.8) +
  geom_rug(data = subset(iber_chr6_res, outa2 == "Candidate SNP"), sides="b") +
  geom_rug(data = subset(iber_chr6_res, outa2 == "Top 10 hits"), sides="b") +
  labs(title="Iberian cline: candidate SNPs for strong selection on homozygotes",
       subtitle="Chr. 6 selection scan conducted on 10 365 SNPs in 109 modern DNA and 393 aDNA samples") +
  xlab("Physical position (in BP units)") + ylab(expression(alpha[2]*' coefficient')) +
  scale_color_manual(values=c("#47abd8","#D3D3D3","#D01B1B")) + 
  theme_bw() + theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) +
  geom_label_repel(data=iber_chr6_res[iber_chr6_res$outa2 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, 
                   xlim  = c(0, 19677984, 33485635, Inf), family = "Bahnschrift", show.legend=FALSE) + 
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  labs(color="SNP type") + 
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold")) + ylim(-50,1600)


table(iber_candidates_a2$sel)
table(iber_top10_a2$sel)
table(iber_top10_a1$sel)


iber_plot_a1_a2 = iber_chr6_res
iber_plot_a1_a2$out = "Other SNPs"
iber_plot_a1_a2$out[iber_plot_a1_a2$SNP %in% iber_top10_a2$SNP] = "Top 10 directional"
iber_plot_a1_a2$out[iber_plot_a1_a2$SNP %in% iber_top10_a1$SNP] = "Top 10 overdominance"


# Plot alpha1 vs alpha2 coefficients for the results
ggplot(data=iber_plot_a1_a2, aes(x=alpha1_mean, y=alpha2_mean, color=out)) + 
  geom_point(alpha=0.15) + geom_point(data = subset(iber_plot_a1_a2, out == "Top 10 directional", alpha=0.3)) + 
  geom_point(data = subset(iber_plot_a1_a2, out == "Top 10 overdominance", alpha=0.3)) + 
  geom_smooth(method=lm , color="#404040", se=FALSE, linetype="dashed") + 
  xlab(expression(alpha[1]*" coefficient")) + ylab(expression(alpha[2]*" coefficient")) +
  labs(color = "SNP type") + scale_color_manual(values=c("#D3D3D3","#47abd8","#D01B1B")) +
  labs(title="Iberian cline: selection coefficients of the scanned SNPs") +
  labs(subtitle="n = 10 365") +
  geom_label_repel(data=iber_plot_a1_a2[iber_plot_a1_a2$outa1 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, show.legend=FALSE) +
  geom_label_repel(data=iber_plot_a1_a2[iber_plot_a1_a2$outa2 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, show.legend=FALSE) +
  theme_bw() + theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"),
        axis.line = element_line(colour = "black"), 
        panel.border = element_blank()) 

# Plot SDs for top hits
a = iber_top10_a1 %>%rename(diff = diff_a1)
a = subset(a, select = -c(alpha2_mean, alpha2_std))

b = subset(iber_top10_a2, select = -c(alpha1_mean, alpha1_std))
b = b %>%rename(alpha1_mean = alpha2_mean)
b = b %>%rename(alpha1_std = alpha2_std)
iber_top10 = rbind(a, b)

ggplot(iber_top10, aes(x=reorder(SNP, -alpha1_mean),alpha1_mean, fill=sel, color=sel))+
  geom_bar(stat="identity", position=position_dodge(), alpha = 0.7) +
  geom_errorbar(aes(ymin=alpha1_mean-alpha1_std, ymax=alpha1_mean+alpha1_std), 
                width=.2, position=position_dodge(.9), alpha=1, color="#404040") + 
  ylab("Respective selection coefficient value") + xlab("SNP") + 
  labs(title="Iberian cline: selection coefficients and standard deviations for the top hits") +
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  scale_color_manual(values=c("#47abd8","#D01B1B")) +
  scale_fill_manual(name="Selection type",labels=c(expression("Direct ("*alpha[2]*" is bigger)"), 
                                                expression(" Overdominant ("*alpha[1]*"is bigger)")), values=c("#47abd8","#D01B1B")) +
  theme_bw() + theme(text=element_text(family="Bahnschrift", size=11), 
    plot.subtitle = element_text(color = "#404040"),
    plot.title = element_text(size = 14, face = "bold", margin=margin(0,0,15,0)),
    axis.line = element_line(colour = "black"), 
    panel.border = element_blank()) + coord_cartesian(ylim=c(0,5500)) + 
    labs(fill='Selection type') 



# Check how many are incommon
steppe_candidates_a1$iber = steppe_candidates_a1$SNP %in% iber_candidates_a1$SNP
steppe_candidates_a2$iber = steppe_candidates_a2$SNP %in% iber_candidates_a2$SNP

table(steppe_candidates_a1$iber)
table(steppe_candidates_a2$iber)

################### Plot frequencies for top hits
setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/Inputs/pruned_inputs")

iber_pruned = read.delim("iber_pruned_input.csv", header = FALSE, sep=",")

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
iber0 = read.delim("iberian_cline0.frq", header=TRUE, sep="")
iber500 = read.delim("iberian_cline500.frq", header=TRUE, sep="")
iber1000 = read.delim("iberian_cline1000.frq", header=TRUE, sep="")
iber1500 = read.delim("iberian_cline1500.frq", header=TRUE, sep="")
iber2000 = read.delim("iberian_cline2000.frq", header=TRUE, sep="")
iber2500 = read.delim("iberian_cline2500.frq", header=TRUE, sep="")
iber3000 = read.delim("iberian_cline3000.frq", header=TRUE, sep="")
iber3500 = read.delim("iberian_cline3500.frq", header=TRUE, sep="")
iber4000 = read.delim("iberian_cline4000.frq", header=TRUE, sep="")
iber4500 = read.delim("iberian_cline4500.frq", header=TRUE, sep="")
iber5000 = read.delim("iberian_cline5000.frq", header=TRUE, sep="")
iber5500 = read.delim("iberian_cline5500.frq", header=TRUE, sep="")
iber6000 = read.delim("iberian_cline6000.frq", header=TRUE, sep="")
iber6500 = read.delim("iberian_cline6500.frq", header=TRUE, sep="")


names(iber0)[5]<-paste("Present")
names(iber500)[5]<-paste("MAF 1 - 499 BP")
names(iber1000)[5]<-paste("MAF 500 - 999 BP")
names(iber1500)[5]<-paste("MAF 1000 - 1499 BP")
names(iber2000)[5]<-paste("MAF 1500 - 1999 BP")
names(iber2500)[5]<-paste("MAF 2000 - 2499 BP")
names(iber3000)[5]<-paste("MAF 2500 - 2999 BP")
names(iber3500)[5]<-paste("MAF 3000 - 3499 BP")
names(iber4000)[5]<-paste("MAF 3500 - 3999 BP")
names(iber4500)[5]<-paste("MAF 4000 - 4499 BP")
names(iber5000)[5]<-paste("MAF 4500 - 4999 BP")
names(iber5500)[5]<-paste("MAF 5000 - 5499 BP")
names(iber6000)[5]<-paste("MAF 5500 - 5999 BP")
names(iber6500)[5]<-paste("MAF 6000 - 6499 BP")


iber_freq = data.frame(iber500$CHR, iber500$SNP, iber500$A1, iber500$A2,
                       iber0$"Present",
                       iber500$"MAF 1 - 499 BP",
                       iber1000$"MAF 500 - 999 BP",
                       iber1500$"MAF 1000 - 1499 BP",
                       iber2000$"MAF 1500 - 1999 BP",
                       iber2500$"MAF 2000 - 2499 BP",
                       iber3000$"MAF 2500 - 2999 BP",
                       iber3500$"MAF 3000 - 3499 BP",
                       iber4000$"MAF 3500 - 3999 BP",
                       iber4500$"MAF 4000 - 4499 BP",
                       iber5000$"MAF 4500 - 4999 BP",
                       iber5500$"MAF 5000 - 5499 BP",
                       iber6000$"MAF 5500 - 5999 BP",
                       iber6500$"MAF 6000 - 6499 BP")


colnames(iber_freq) <- c("CHR", "SNP", "A1", "A2",
                         "Present",
                         "MAF 1 - 499 BP",
                         "MAF 500 - 999 BP",
                         "MAF 1000 - 1499 BP",
                         "MAF 1500 - 1999 BP",
                         "MAF 2000 - 2499 BP",
                         "MAF 2500 - 2999 BP",
                         "MAF 3000 - 3499 BP",
                         "MAF 3500 - 3999 BP",
                         "MAF 4000 - 4499 BP",
                         "MAF 4500 - 4999 BP",
                         "MAF 5000 - 5499 BP",
                         "MAF 5500 - 5999 BP",
                         "MAF 6000 - 6499 BP")

iber_chr6_freq = iber_freq[iber_freq$SNP %in% iber_top10_a1$SNP,]
iber_chr6_freq = iber_chr6_freq[,c(2,5:18)]
iber_chr6_freq = melt(iber_chr6_freq ,  id.vars = 'SNP', variable.name = 'Period')


ggplot(iber_chr6_freq, aes(Period, value, group=SNP)) +
  geom_line(size=0.01, aes(color=SNP)) + 
  scale_x_discrete(limits = rev(levels(iber_chr6_freq$Period))) +
  theme_bw()

