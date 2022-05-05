### 5. Generating frequency bins for the analysis ###

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/Inputs/pruned_inputs")
steppe_input = read.table("steppe_pruned_input.csv", header = TRUE, sep=",")

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
steppe_modern_freq = read.table("steppe_cline0.frq", header = TRUE)

steppe_modern_freq = steppe_modern_freq[ , c(1,2,5)]
steppe_modern_freq = steppe_modern_freq[complete.cases(steppe_modern_freq), ]
steppe_modern_freq = steppe_modern_freq[steppe_modern_freq$SNP %in% steppe_input$SNP,]

steppe_modern_freq =  steppe_modern_freq %>% 
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
steppe_neutrality_snps = steppe_modern_freq %>% group_by(freq_bin) %>% slice_sample(n=100)
steppe_neutrality_input = steppe_input[steppe_input$SNP %in% steppe_neutrality_snps$SNP, ]

setwd("C:/Users/ACER/Desktop")
fwrite(steppe_neutrality_input, "steppe_neutrality_input.csv")


#Read the results for bins' analysis
steppe_neutrality_results = read.table("steppe_neutrality_results.csv", header = FALSE, sep=",")
colnames(steppe_neutrality_results) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")

steppe_neutrality_snps = steppe_modern_freq[,c(1,2,4)]
steppe_neutrality_results = merge(steppe_neutrality_results, steppe_modern_freq, by="SNP")

#Check and plot the distributions of alpha1 and alpha2 
setDT(steppe_neutrality_results)
steppe_neutr_a1 = as.data.frame(steppe_neutrality_results[, as.list(summary(alpha1_mean)), by = freq_bin])
steppe_neutr_a2 = as.data.frame(steppe_neutrality_results[, as.list(summary(alpha2_mean)), by = freq_bin])


steppe_neutr_a1$freq_bin = factor(steppe_neutr_a1$freq_bin,
                                  levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                           "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                           "85-90%","90-95%","95-100%"))

steppe_neutr_a2$freq_bin = factor(steppe_neutr_a2$freq_bin,
                                  levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                           "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                           "85-90%","90-95%","95-100%"))



#alpha1 mean in different bins barplot
ggplot(data=steppe_neutr_a1, aes(x=freq_bin, y=Median)) + geom_bar(stat="identity")

#alpha2 mean in different bins barplot
ggplot(data=steppe_neutr_a2, aes(x=freq_bin, y=Median)) + geom_bar(stat="identity")

#Plot distribution of alpha1 and alpha2 within the bins
steppe_neutr_a1$Coefficient = "Alpha 1" 
steppe_neutr_a2$Coefficient = "Alpha 2"
steppe_neutr_merged = rbindlist(list(steppe_neutr_a1,steppe_neutr_a2))

ggplot(data = steppe_neutr_merged, aes(x = freq_bin, y = Median, color = Coefficient, fill = Coefficient)) + 
  geom_bar(stat = "identity", alpha=0.5) + scale_color_manual(values = c("#E7B800","#00AFBB")) +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) + xlab("Frequency bin") + 
  scale_x_discrete(breaks = levels(steppe_neutr_merged$freq_bin)[c(T, rep(F, 8))])


#Plot alpha1 vs alpha2
setwd("C:/Users/ACER/Desktop")
steppe_neutrality_results = read.table("steppe_neutrality_results.csv", header = FALSE, sep=",")
colnames(steppe_neutrality_results) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")
steppe_neutrality_results = steppe_neutrality_results[-c(1403,137,133,851,827,311,31, 1057,415), ]

steppe_neutrality_results$Selection[steppe_neutrality_results$alpha2_mean > steppe_neutrality_results$alpha1_mean & steppe_neutrality_results$alpha1_mean > 0] = "Directional"
steppe_neutrality_results$Selection[steppe_neutrality_results$alpha1_mean > steppe_neutrality_results$alpha2_mean & steppe_neutrality_results$alpha1_mean > 0] = "Overdominant"
steppe_neutrality_results$Selection[steppe_neutrality_results$alpha2_mean > steppe_neutrality_results$alpha1_mean & 0 > steppe_neutrality_results$alpha1_mean] = "Underdominant"

table(steppe_neutrality_results$Selection)

ggplot(data=steppe_neutrality_results, aes(x=alpha1_mean, y=alpha2_mean, color=Selection)) + 
  geom_point() + geom_smooth(method=lm , color="black", se=FALSE) + 
  xlab(expression(alpha[1]*" coefficient")) + ylab(expression(alpha[2]*" coefficient")) 


#Read the results for the entire chromosome 6
setwd("C:/Users/ACER/Desktop")
steppe_chr6_res = read.table("steppe_chr6._results.csv", header = FALSE, sep=",")
colnames(steppe_chr6_res) = c("SNP","alpha1_mean","alpha1_std", "alpha2_mean", "alpha2_std")



setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
snp_position = read.table("steppe_cline0.map")
snp_position = snp_position[, c(1:4)]
colnames(snp_position) = c("CHR","SNP","gpos", "phpos")
steppe_chr6_snp = snp_position[snp_position$SNP %in% steppe_chr6_res$SNP, ]
steppe_chr6_res = merge(steppe_chr6_snp, steppe_chr6_res, by="SNP")

steppe_modern_freq = steppe_modern_freq[,c(2,4)]
steppe_chr6_res = merge(steppe_chr6_res, steppe_modern_freq, by="SNP")

steppe_chr6_res$freq_bin = factor(steppe_chr6_res$freq_bin,
                                  levels=c("0-5%","5-10%","10-15%","15-20%","20-25%","25-30%","30-35%","35-40%","40-45%",
                                           "45-50%","50-55%","55-60%","60-65%","65-70%","70-75%","75-80%","80-85%",
                                           "85-90%","90-95%","95-100%"))

#Look for outliers in the entire data set

#Calculate IQR
steppe_neutr_a1$IQR = steppe_neutr_a1$`3rd Qu.` - steppe_neutr_a1$`1st Qu.`
steppe_neutr_a2$IQR = steppe_neutr_a2$`3rd Qu.` - steppe_neutr_a2$`1st Qu.`

#Calculate upper fence for strong outliers
steppe_neutr_a1$UFSO = steppe_neutr_a1$`3rd Qu.` + 3 * steppe_neutr_a1$IQR
steppe_neutr_a2$UFSO = steppe_neutr_a2$`3rd Qu.` + 3 * steppe_neutr_a2$IQR

#Calculate lower fence for strong outliers
steppe_neutr_a1$LFSO = steppe_neutr_a1$`1st Qu.` - 3 * steppe_neutr_a1$IQR
steppe_neutr_a2$LFSO = steppe_neutr_a2$`1st Qu.` - 3 * steppe_neutr_a2$IQR

#UF for a1

steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="0-5%"] = 38.00682 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="5-10%"] = 66.77780 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="10-15%"] = 57.26574 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="15-20%"] = 60.35344 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="20-25%"] = 75.20785 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="25-30%"] = 85.31126 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="30-35%"] = 103.99833 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="35-40%"] = 99.65047 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="40-45%"] = 109.39255 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="45-50%"] = 146.00248 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="50-55%"] = 92.59801 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="55-60%"] = 207.85205 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="60-65%"] = 617.29670 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="65-70%"] = 1528.41373 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="70-75%"] = 1670.60977 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="75-80%"] = 1778.59855 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="80-85%"] = 1694.70605 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="85-90%"] = 1234.77621 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="90-95%"] = 1486.04851 
steppe_chr6_res$UFSO.a1[steppe_chr6_res$freq_bin=="95-100%"] = 1554.31782 

#UF for a2
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="0-5%"] = 36.92502
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="5-10%"] = 63.86843
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="10-15%"] = 41.58715
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="15-20%"] = 49.23797
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="20-25%"] = 39.60533
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="25-30%"] = 37.69548
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="30-35%"] = 39.23020
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="35-40%"] = 32.71402
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="40-45%"] = 41.07724
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="45-50%"] = 75.15424
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="50-55%"] = 101.06842
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="55-60%"] = 234.77746
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="60-65%"] = 674.96961
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="65-70%"] = 1448.52875
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="70-75%"] = 1444.64570
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="75-80%"] = 1547.44540
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="80-85%"] = 1532.99080
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="85-90%"] = 1119.19094
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="90-95%"] = 1442.28504
steppe_chr6_res$UFSO.a2[steppe_chr6_res$freq_bin=="95-100%"] = 	1355.40018

#Find & classify a1 outliers

steppe_chr6_res$outa1[steppe_chr6_res$UFSO.a1 > steppe_chr6_res$alpha1_mean | steppe_chr6_res$alpha1_mean < steppe_chr6_res$alpha2_mean] = "Regular SNP"
steppe_chr6_res$outa1[steppe_chr6_res$UFSO.a1 < steppe_chr6_res$alpha1_mean & steppe_chr6_res$alpha1_mean > steppe_chr6_res$alpha2_mean] = "Candidate SNP"


table(steppe_chr6_res$outa1)

d = subset(steppe_chr6_res, outa1=="Candidate SNP" & alpha1_mean > alpha2_mean)

#Find & classify a2 outliers

steppe_chr6_res$outa2[steppe_chr6_res$UFSO.a2 > steppe_chr6_res$alpha2_mean | steppe_chr6_res$alpha1_mean > steppe_chr6_res$alpha2_mean] = "Regular SNP"
steppe_chr6_res$outa2[steppe_chr6_res$UFSO.a2 < steppe_chr6_res$alpha2_mean & steppe_chr6_res$alpha1_mean < steppe_chr6_res$alpha2_mean] = "Candidate SNP"

table(steppe_chr6_res$outa2)

e = subset(steppe_chr6_res, outa2=="Candidate SNP" & alpha1_mean < alpha2_mean)

#Estimate the amount of selection candidates within and outside HLA (for Chi-square tests)

hla_out_a1 = subset(d, phpos >= 29677984 & 33485635 >= phpos)
hla_non_out_a1 = subset(d, phpos < 29677984 | 33485635 < phpos)

hla_out_a2 = subset(e, phpos >= 29677984 & 33485635 >= phpos)
hla_non_out_a2 = subset(e, phpos < 29677984 | 33485635 < phpos)

hla_SNPs = subset(steppe_chr6_res, phpos >= 29677984 & 33485635 >= phpos)
non_hla_SNPs = subset(steppe_chr6_res, phpos < 29677984 | 33485635 < phpos)

#Classify selection type

steppe_chr6_res$sel[steppe_chr6_res$alpha2_mean > steppe_chr6_res$alpha1_mean & steppe_chr6_res$alpha1_mean > 0] = "Direct"
steppe_chr6_res$sel[steppe_chr6_res$alpha1_mean > steppe_chr6_res$alpha2_mean & steppe_chr6_res$alpha1_mean > 0] = "Overdominant"
steppe_chr6_res$sel[steppe_chr6_res$alpha2_mean > steppe_chr6_res$alpha1_mean & 0 > steppe_chr6_res$alpha1_mean] = "Underdominant"

table(steppe_chr6_res$sel)



# Find top hits for alpha1
steppe_candidates_a1 = subset(steppe_chr6_res, outa1=="Candidate SNP" & alpha1_mean > alpha2_mean)
steppe_candidates_a1$diff_a1 = steppe_candidates_a1$alpha1_mean/steppe_candidates_a1$UFSO.a1
steppe_candidates_a1 = steppe_candidates_a1[order(steppe_candidates_a1$diff_a1, decreasing = TRUE), ] 

steppe_top10_a1 = top_n(steppe_candidates_a1, 10, diff_a1)
steppe_chr6_res$outa1[steppe_chr6_res$SNP %in% steppe_top10_a1$SNP] = "Top 10 hits"

#Export the results
#steppe_top10_export = steppe_top10_a1[,c(1,4,5,14)]
#fwrite(steppe_top10_export, "steppe_top10.csv")

hla_steppe_a1 = subset(steppe_candidates_a1, phpos >= 29677984 & 33485635 >= phpos)
non_hla_steppe_a1 = subset(steppe_candidates_a1, phpos < 29677984 | 33485635 < phpos)

ggplot(steppe_chr6_res, aes(x=phpos, y=alpha1_mean, color=outa1)) + 
  geom_point(alpha=0.15) +
  geom_vline(xintercept = 29677984, alpha=0.2, linetype="longdash") +
  geom_vline(xintercept = 33485635, alpha=0.2, linetype="longdash") +
  geom_point(data = subset(steppe_chr6_res, outa1 == "Candidate SNP"), alpha=0.6) + 
  geom_point(data = subset(steppe_chr6_res, outa1 == "Top 10 hits"), alpha=0.8) +
  geom_rug(data = subset(steppe_chr6_res, outa1 == "Candidate SNP"), sides="b") +
  geom_rug(data = subset(steppe_chr6_res, outa1 == "Top 10 hits"), sides="b") +
  labs(title="Steppe cline: candidate SNPs for strong selection on heterozygotes",
       subtitle="Chr. 6 selection scan conducted on 8381 SNPs in 126 modern DNA and 727 aDNA samples") +
  xlab("Physical position (in BP units)") + ylab(expression(alpha[1]*' coefficient')) +
  scale_color_manual(values=c("#ffa600", "#D3D3D3", "#bc5090")) + 
   theme_bw() + theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) +
  geom_label_repel(data=steppe_chr6_res[steppe_chr6_res$outa1 == "Top 10 hits",],
      aes(label=SNP), size=2.5, min.segment.length = 0, 
      xlim  = c(0, 19677984, 33485635, Inf), family = "Bahnschrift", show.legend=FALSE) + 
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  labs(color="SNP type") + 
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"))


table(steppe_candidates_a1$sel)


# Find top hits for alpha2
steppe_candidates_a2 = subset(steppe_chr6_res, outa2=="Candidate SNP" & alpha1_mean < alpha2_mean)
steppe_candidates_a2$diff = steppe_candidates_a2$alpha2_mean/steppe_candidates_a2$UFSO.a2
steppe_candidates_a2 = steppe_candidates_a2[order(steppe_candidates_a2$diff, decreasing = TRUE), ] 

steppe_top10_a2 = top_n(steppe_candidates_a2, 10, diff)
steppe_chr6_res$outa2[steppe_chr6_res$SNP %in% steppe_top10_a2$SNP] = "Top 10 hits"

hla_steppe_a2 = subset(steppe_candidates_a2, phpos >= 29677984 & 33485635 >= phpos)

ggplot(steppe_chr6_res, aes(x=phpos, y=alpha2_mean, color=outa2)) + 
  geom_point(alpha=0.15) +
  geom_vline(xintercept = 29677984, alpha=0.2, linetype="longdash") +
  geom_vline(xintercept = 33485635, alpha=0.2, linetype="longdash") +
  geom_point(data = subset(steppe_chr6_res, outa2 == "Candidate SNP"), alpha=0.6) + 
  geom_point(data = subset(steppe_chr6_res, outa2 == "Top 10 hits"), alpha=0.8) +
  geom_rug(data = subset(steppe_chr6_res, outa2 == "Candidate SNP"), sides="b") +
  geom_rug(data = subset(steppe_chr6_res, outa2 == "Top 10 hits"), sides="b") +
  labs(title="Steppe cline: candidate SNPs for strong selection on homozygotes",
       subtitle="Chr. 6 selection scan conducted on 8381 SNPs in 126 modern DNA and 727 aDNA samples") +
  xlab("Physical position (in BP units)") + ylab(expression(alpha[2]*' coefficient')) +
  scale_color_manual(values=c("#ffa600", "#D3D3D3", "#bc5090")) + 
  theme_bw() + theme(axis.line = element_line(colour = "black"), panel.border = element_blank()) +
  geom_label_repel(data=steppe_chr6_res[steppe_chr6_res$outa2 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, 
                   xlim  = c(0, 19677984, 33485635, Inf), family = "Bahnschrift", show.legend=FALSE) + 
  guides(fill = guide_legend(override.aes = aes(label = ""))) +
  labs(color="SNP type") + ylim(-50,1700) +
  theme(text=element_text(family="Bahnschrift", size=11), 
        plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
        plot.title = element_text(size = 14, face = "bold"))

table(steppe_candidates_a2$sel)


table(steppe_top10_a2$sel)
table(steppe_top10_a1$sel)


steppe_plot_a1_a2 = steppe_chr6_res
steppe_plot_a1_a2$out = "Other SNPs"
steppe_plot_a1_a2$out[steppe_plot_a1_a2$SNP %in% steppe_top10_a2$SNP] = "Top 10 directional"
steppe_plot_a1_a2$out[steppe_plot_a1_a2$SNP %in% steppe_top10_a1$SNP] = "Top 10 overdominance"


# Plot alpha1 vs alpha2 coefficients for the results
ggplot(data=steppe_plot_a1_a2, aes(x=alpha1_mean, y=alpha2_mean, color=out)) + 
  geom_point(alpha=0.15) + geom_point(data = subset(steppe_plot_a1_a2, out == "Top 10 directional", alpha=0.3)) + 
  geom_point(data = subset(steppe_plot_a1_a2, out == "Top 10 overdominance", alpha=0.3)) + 
  geom_smooth(method=lm , color="#404040", se=FALSE, linetype="dashed") + 
  xlab(expression(alpha[1]*" coefficient")) + ylab(expression(alpha[2]*" coefficient")) +
  labs(color = "SNP type") + scale_color_manual(values=c("#D3D3D3","#ffa600","#bc5090")) +
  labs(title="Steppe cline: selection coefficients of the scanned SNPs") +
  labs(subtitle="n = 8 381") +
  geom_label_repel(data=steppe_plot_a1_a2[steppe_plot_a1_a2$outa1 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, show.legend=FALSE) +
  geom_label_repel(data=steppe_plot_a1_a2[steppe_plot_a1_a2$outa2 == "Top 10 hits",],
                   aes(label=SNP), size=2.5, min.segment.length = 0, show.legend=FALSE) +
  theme_bw() + theme(text=element_text(family="Bahnschrift", size=11), 
                     plot.subtitle = element_text(color = "#404040", margin=margin(0,0,15,0)),
                     plot.title = element_text(size = 14, face = "bold"),
                     axis.line = element_line(colour = "black"), 
                     panel.border = element_blank()) 

ggplot(data=steppe_plot_a1_a2, aes(x=alpha1_mean, y=alpha2_mean, color=out)) + 
  geom_point(alpha=0.15) + 
  geom_point(data = subset(steppe_plot_a1_a2, out == "Top 10 directional", alpha=1)) + 
  geom_point(data = subset(steppe_plot_a1_a2, out == "Top 10 overdominance", alpha=1)) + 
  geom_smooth(method=lm , color="black", se=FALSE) + 
  xlab(expression(alpha[1]*" coefficient")) + ylab(expression(alpha[2]*" coefficient"))+
  labs(color = "SNP type") + scale_color_manual(values=c("#D3D3D3","#ffca00", "#8B008B"))+
  theme_bw() + labs(title=expression("Steppe cline "*alpha[1]*" vs "*alpha[2]))

# Plot SDs for top hits
a = steppe_top10_a1 %>%rename(diff = diff_a1)
a = subset(a, select = -c(alpha2_mean, alpha2_std))

b = subset(steppe_top10_a2, select = -c(alpha1_mean, alpha1_std))
b = b %>%rename(alpha1_mean = alpha2_mean)
b = b %>%rename(alpha1_std = alpha2_std)
steppe_top10 = rbind(a, b)


ggplot(steppe_top10, aes(x=reorder(SNP, -alpha1_mean),alpha1_mean, fill=sel, color=sel))+
  geom_bar(stat="identity", position=position_dodge(), alpha = 0.8) +
  geom_errorbar(aes(ymin=alpha1_mean-alpha1_std, ymax=alpha1_mean+alpha1_std), 
                width=.2, position=position_dodge(.9), alpha=1, color="#404040") + 
  ylab("Respective selection coefficient value") + xlab("SNP") + 
  labs(title="Steppe cline: selection coefficients and standard deviations for the top hits") +
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  scale_color_manual(values=c("#ffca00", "#8B008B")) +
  scale_fill_manual(name="Selection type",
                    labels=c(expression("Direct ("*alpha[2]*" is bigger)"), 
                            expression(" Overdominant ("*alpha[1]*"is bigger)")), 
                    values=c("#ffa600","#bc5090")) +
  theme_bw() + theme(text=element_text(family="Bahnschrift", size=11), 
                     plot.subtitle = element_text(color = "#404040"),
                     plot.title = element_text(size = 14, face = "bold", margin=margin(0,0,15,0)),
                     axis.line = element_line(colour = "black"), 
                     panel.border = element_blank()) + coord_cartesian(ylim=c(0,4000)) + 
  labs(fill='Selection type') 



################ Plot frequency 
setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/Inputs/pruned_inputs")

steppe_pruned = read.delim("steppe_pruned_input.csv", header = FALSE, sep=",")

setwd("C:/Users/ACER/Desktop/Undergraduate dissertation/v44.3_1240K_public")
steppe0 = read.delim("steppe_cline0.frq", header=TRUE, sep="")
steppe500 = read.delim("steppe_cline500.frq", header=TRUE, sep="")
steppe1000 = read.delim("steppe_cline1000.frq", header=TRUE, sep="")
steppe1500 = read.delim("steppe_cline1500.frq", header=TRUE, sep="")
steppe2000 = read.delim("steppe_cline2000.frq", header=TRUE, sep="")
steppe2500 = read.delim("steppe_cline2500.frq", header=TRUE, sep="")
steppe3000 = read.delim("steppe_cline3000.frq", header=TRUE, sep="")
steppe3500 = read.delim("steppe_cline3500.frq", header=TRUE, sep="")
steppe4000 = read.delim("steppe_cline4000.frq", header=TRUE, sep="")
steppe4500 = read.delim("steppe_cline4500.frq", header=TRUE, sep="")
steppe5000 = read.delim("steppe_cline5000.frq", header=TRUE, sep="")
steppe5500 = read.delim("steppe_cline5500.frq", header=TRUE, sep="")
steppe6000 = read.delim("steppe_cline6000.frq", header=TRUE, sep="")
steppe6500 = read.delim("steppe_cline6500.frq", header=TRUE, sep="")


names(steppe0)[5]<-paste("Present")
names(steppe500)[5]<-paste("MAF 1 - 499 BP")
names(steppe1000)[5]<-paste("MAF 500 - 999 BP")
names(steppe1500)[5]<-paste("MAF 1000 - 1499 BP")
names(steppe2000)[5]<-paste("MAF 1500 - 1999 BP")
names(steppe2500)[5]<-paste("MAF 2000 - 2499 BP")
names(steppe3000)[5]<-paste("MAF 2500 - 2999 BP")
names(steppe3500)[5]<-paste("MAF 3000 - 3499 BP")
names(steppe4000)[5]<-paste("MAF 3500 - 3999 BP")
names(steppe4500)[5]<-paste("MAF 4000 - 4499 BP")
names(steppe5000)[5]<-paste("MAF 4500 - 4999 BP")
names(steppe5500)[5]<-paste("MAF 5000 - 5499 BP")
names(steppe6000)[5]<-paste("MAF 5500 - 5999 BP")
names(steppe6500)[5]<-paste("MAF 6000 - 6499 BP")


steppe_freq = data.frame(steppe500$CHR, steppe500$SNP, steppe500$A1, steppe500$A2,
                       steppe0$"Present",
                       steppe500$"MAF 1 - 499 BP",
                       steppe1000$"MAF 500 - 999 BP",
                       steppe1500$"MAF 1000 - 1499 BP",
                       steppe2000$"MAF 1500 - 1999 BP",
                       steppe2500$"MAF 2000 - 2499 BP",
                       steppe3000$"MAF 2500 - 2999 BP",
                       steppe3500$"MAF 3000 - 3499 BP",
                       steppe4000$"MAF 3500 - 3999 BP",
                       steppe4500$"MAF 4000 - 4499 BP",
                       steppe5000$"MAF 4500 - 4999 BP",
                       steppe5500$"MAF 5000 - 5499 BP",
                       steppe6000$"MAF 5500 - 5999 BP",
                       steppe6500$"MAF 6000 - 6499 BP")


colnames(steppe_freq) <- c("CHR", "SNP", "A1", "A2",
                         "Present",
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
                         "6000 - 6499 BP")

steppe_best_SNPs = subset(steppe_top10, SNP=="rs7756516" | SNP=="rs3130062" | SNP=="rs3132453")

steppe_chr6_freq = steppe_freq[steppe_freq$SNP %in% steppe_best_SNPs$SNP,]
steppe_chr6_freq = steppe_chr6_freq[,c(2,5:18)]
steppe_chr6_freq = melt(steppe_chr6_freq ,  id.vars = 'SNP', variable.name = 'Period')


ggplot(steppe_chr6_freq, aes(Period, value, group=SNP)) +
  geom_line(size=0.01, aes(color=SNP)) + 
  scale_x_discrete(limits = rev(levels(steppe_chr6_freq$Period))) +
  theme_bw()


steppe_best_SNPs = subset(steppe_top10, SNP=="rs7756516")
steppe_chr6_freq = steppe_freq[steppe_freq$SNP %in% steppe_best_SNPs$SNP,]
steppe_chr6_freq = steppe_chr6_freq[,c(2,5:18)]
steppe_chr6_freq = melt(steppe_chr6_freq ,  id.vars = 'SNP', variable.name = 'Period')

ggplot(steppe_chr6_freq, aes(Period, value, group=SNP)) +
  geom_point() + ylab("Frequency") +
  scale_x_discrete(limits = rev(levels(steppe_chr6_freq$Period)), guide = guide_axis(angle = 45)) +
   stat_smooth(geom='ribbon', aes(ymin = ifelse(..ymin.. < 0, 0, ..ymin..)), alpha = .3) +
   geom_smooth(se = FALSE) + theme_bw() + theme_1()

