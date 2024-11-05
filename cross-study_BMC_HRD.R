library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(zoo)
library(scales)
library(data.table)
library(ggrepel)
library(forcats)
library(readxl)
library(scales)
library(ggridges)
library(survival)
library(wesanderson)
library(viridis)
library(pals)

BRCA<-c("BRCA1", "BRCA2")

HRR_13<-c( "ATM", "BRIP1", "PALB2", "RAD51C", "BARD1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PPP2R2A", "RAD51B", "RAD51D", "RAD54L")

HRR_myr<-c("BRCA1", "BRCA2", "ATM", "BRIP1", "PALB2", "RAD51C", "BARD1", "CDK12", "CHEK1", "CHEK2", "FANCL", "PPP2R2A", "RAD51B", "RAD51D", "RAD54L")


###read in all files:
PAOLA<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/PAOLA/PAOLA_summary_paper_withfails.csv")
OPINION<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/OPINION/OPINION_paper_refined.csv")
LIGHT<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/LIGHT/LIGHT_paper_refined.csv")
Study19<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Study19/AZO_HRD_HRR_results_summary_Ecode_added_nonHRR_GS_added.csv")
SOLO1<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/SOLO1/SOLO1_Myriad_data_summary_paper.csv")
SOLO2<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/SOLO2/SOLO2_Myriad_data_summary_20200624.csv")



####remove test fails:
OPINION<-OPINION[!OPINION$Status =="Test Failed/Cancelled",]
LIGHT<-LIGHT[!LIGHT$Status =="Test Failed/Cancelled",]




#HRD fails
length(PAOLA$SUBJECT[PAOLA$HRD_Score =="FAILED"])
length(PAOLA$SUBJECT[PAOLA$HRD_Score =="FAILED"])/length(PAOLA$SUBJECT)*100
PAOLA_fail<-c(length(PAOLA$SUBJECT[PAOLA$HRD_Score =="FAILED"]),length(PAOLA$SUBJECT[PAOLA$HRD_Score =="FAILED"])/length(PAOLA$SUBJECT)*100)

length(OPINION$SUBJECT[OPINION$HRD_Score =="FAILED"])
length(OPINION$SUBJECT[OPINION$HRD_Score =="FAILED"])/length(OPINION$SUBJECT)
OPINION_fail<-c(length(OPINION$SUBJECT[OPINION$HRD_Score =="FAILED"]),length(OPINION$SUBJECT[OPINION$HRD_Score =="FAILED"])/length(OPINION$SUBJECT)*100)


length(LIGHT$SUBJECT[LIGHT$HRD_Score =="FAILED"])
length(LIGHT$SUBJECT[LIGHT$HRD_Score =="FAILED"])/length(LIGHT$SUBJECT)*100
LIGHT_fail<-c(length(LIGHT$SUBJECT[LIGHT$HRD_Score =="FAILED"]),length(LIGHT$SUBJECT[LIGHT$HRD_Score =="FAILED"])/length(LIGHT$SUBJECT)*100)


length(Study19$SUBJECT[Study19$HRD_Score =="FAILED"])
length(Study19$SUBJECT[Study19$HRD_Score =="FAILED"])/length(Study19$SUBJECT)*100

Study19_fail<-c(length(Study19$SUBJECT[Study19$HRD_Score =="FAILED"]),length(Study19$SUBJECT[Study19$HRD_Score =="FAILED"])/length(Study19$SUBJECT)*100)


length(SOLO1$SUBJECT[SOLO1$HRD_Score =="FAILED"])
length(SOLO1$SUBJECT[SOLO1$HRD_Score =="FAILED"])/length(SOLO1$SUBJECT)*100

SOLO1_fail<-c(length(SOLO1$SUBJECT[SOLO1$HRD_Score =="FAILED"]),length(SOLO1$SUBJECT[SOLO1$HRD_Score =="FAILED"])/length(SOLO1$SUBJECT)*100)


length(SOLO2$SUBJECT[SOLO2$HRD_Score =="FAILED"])
length(SOLO2$SUBJECT[SOLO2$HRD_Score =="FAILED"])/length(SOLO2$SUBJECT)*100

SOLO2_fail<-c(length(SOLO2$SUBJECT[SOLO2$HRD_Score =="FAILED"]),length(SOLO2$SUBJECT[SOLO2$HRD_Score =="FAILED"])/length(SOLO2$SUBJECT)*100)

all_HRD_fail<-rbind(PAOLA_fail,OPINION_fail,LIGHT_fail,Study19_fail,SOLO1_fail,SOLO2_fail)

all_HRD<-rbind(PAOLA,OPINION,LIGHT,Study19,SOLO1,SOLO2)
length(all_HRD$SUBJECT[all_HRD$HRD_Score =="FAILED"])
length(all_HRD$SUBJECT[all_HRD$HRD_Score !="FAILED"])

length(all_HRD$SUBJECT[all_HRD$HRD_Score =="FAILED"])/length(all_HRD$SUBJECT)*100

all_HRD_fail_reason<-data.frame(table(all_HRD$Fail_Reason[all_HRD$HRD_Score =="FAILED"]))
all_HRD_fail_reason$Percent <-(all_HRD_fail_reason$Freq/sum(all_HRD_fail_reason$Freq))*100

write.csv(all_HRD,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20200619.csv",row.names=FALSE)


all_HRD<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20200619.csv")
all_HRD$Population<-ifelse(all_HRD$Study =="PAOLA" | all_HRD$Study =="SOLO1", "CR/PR after 1L Platinum Based Chemotherapy", "PSR")
all_HRD$Variant_Subtype<-gsub("frameshift/indel","Frameshift/Indel",all_HRD$Variant_Subtype)
all_HRD$Variant_Subtype<-gsub("Loss/Rearrangment","Loss/Rearrangement",all_HRD$Variant_Subtype)

write.csv(all_HRD,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20200629.csv",row.names=FALSE)

all_HRD<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20200629.csv")

all_HRD$Status_3<-ifelse(all_HRD$Status =="tBRCAm", "tBRCAm", "Non-tBRCAm")
all_HRD$HRD_status<-ifelse(all_HRD$Status =="tBRCAm","HRD positive (tBRCAm)",
                           ifelse(all_HRD$HRD_Score >=42 & all_HRD$Status == "HRRm and non-BRCAm","HRD positive (HRRm)",
                                  ifelse(all_HRD$HRD_Score >=42,"HRD positive (non-HRRm)",
                                         ifelse(all_HRD$HRD_Score < 42 & all_HRD$Status == "HRRm and non-BRCAm", "HRD negative (HRRm)",
                                                ifelse(all_HRD$HRD_Score < 42, "HRD negative (non-HRRm)",
                                                       ifelse(all_HRD$HRD_Score =="FAILED" & all_HRD$Status == "HRRm and non-BRCAm", "HRD unknown (HRRm)","HRD unknown"))))))

write.csv(all_HRD,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20201020_int.csv",row.names=FALSE)
all_HRD<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20201020_int.csv")

  
all_HRD_plots<-all_HRD[!all_HRD$HRD_Score =="FAILED",]
all_HRD_plots<-all_HRD_plots[!is.na(all_HRD_plots$Status),]

#all_HRD_plots<-all_HRD_plots[!is.na(all_HRD$HRD_Score),]

all_HRD_plots$HRD_Score <-as.numeric(as.character(all_HRD_plots$HRD_Score))

all_HRD_plots$Status <- gsub("HRRm and non-BRCAm","Non-BRCA HRRm", all_HRD_plots$Status)
  
all_HRD_plots$Status <- factor(all_HRD_plots$Status, levels=c( "tBRCAm","Non-BRCA HRRm","Non-HRRm")) #reordering zygosity variables to appear in this order in the barplot
all_HRD_plots$Status_3 <- factor(all_HRD_plots$Status_3, levels=c( "tBRCAm","Non-tBRCAm")) #reordering zygosity variables to appear in this order in the barplot

all_HRD_plots$Study <- factor(all_HRD_plots$Study, levels=c( "PAOLA","OPINION","LIGHT","Study 19", "SOLO1", "SOLO2")) #reordering zygosity variables to appear in this order in the barplot
all_HRD_plots$key<-paste(all_HRD_plots$SUBJECT,all_HRD_plots$Study,sep="_")

t.test(all_HRD_plots$HRD_Score[all_HRD_plots$Status =="tBRCAm"],all_HRD_plots$HRD_Score[all_HRD_plots$Status =="Non-HRRm"])
t.test(all_HRD_plots$HRD_Score[all_HRD_plots$Status =="tBRCAm"],all_HRD_plots$HRD_Score[all_HRD_plots$Status =="Non-BRCA HRRm"])




all_HRD_plots<-all_HRD_plots[!duplicated(all_HRD_plots$key),]
###Proportion of tBRCAm tumours >42
#all_HRD$HRD_Score <-as.numeric(as.character(all_HRD$HRD_Score))

dim(all_HRD[all_HRD$Status =="tBRCAm",])

dim(all_HRD[all_HRD$Status =="tBRCAm",])/dim(all_HRD)

dim(all_HRD[all_HRD$Status !="tBRCAm",])
dim(all_HRD[all_HRD$Status !="tBRCAm",])/dim(all_HRD)

dim(all_HRD[all_HRD$Status =="tBRCAm" & all_HRD$HRD_Score != "FAILED",])
dim(all_HRD[all_HRD$Status =="tBRCAm" & all_HRD$HRD_Score != "FAILED",])/dim(all_HRD[all_HRD$Status =="tBRCAm",])


dim(all_HRD[all_HRD$Status =="tBRCAm"  & all_HRD$HRD_Score != "FAILED" & as.numeric(as.character(all_HRD$HRD_Score)) >= 42,])

dim(all_HRD[all_HRD$Status =="tBRCAm"  & all_HRD$HRD_Score != "FAILED" & as.numeric(as.character(all_HRD$HRD_Score)) >= 33,])


#####histogram with density plots:
tiff("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Figures/Plot1_density_hist_all_studies.tiff", units="in", width=5, height=5, res=300)
ggplot(all_HRD_plots, aes(x=HRD_Score))+
  geom_histogram(aes(y = ..density..),color="blue",fill="white",binwidth=2)+
  geom_density(alpha=.2,fill="#FF6655")+
  xlab("HRD Score")+
  ylab("Normalised Density")+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 101), breaks = seq(0, 100, by = 25))+
  geom_vline(xintercept=42,linetype="dashed")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8,face="bold"))+
  theme(legend.title=element_text(size=8,face="bold"))+
  theme(legend.text=element_text(size=8))+
  theme(legend.position = 'top')
dev.off()


####density plots: ridges

all_HRD_plots$Status_Ridge<-case_when(all_HRD_plots$Status =="tBRCAm" ~"tBRCAm (n=863)",
                                      all_HRD_plots$Status =="Non-BRCA HRRm" ~"Non-BRCA HRRm (n=107)",
                                      all_HRD_plots$Status =="Non-HRRm" ~"Non-HRRm (n=868)")

all_HRD_plots$Status_Ridge <- factor(all_HRD_plots$Status_Ridge, levels=c( "tBRCAm (n=863)","Non-BRCA HRRm (n=107)","Non-HRRm (n=868)")) #reordering zygosity variables to appear in this order in the barplot

ggplot(all_HRD_plots,aes(y=Status_Ridge))+
  geom_density_ridges(scale= 4,
                      aes(x = HRD_Score, fill =Status_Ridge),
                      alpha = .7, color = NA, from = 0, to=100) +
  labs(
    x = "GIS",
    y = "",
    title = "",
    subtitle = "",
    caption = ""
  )+
  #facet_wrap(~Study,scales="free_x",ncol=1) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, 1.8)))+ 
  #scale_y_discrete(expand = c(0.1, -0.1))+
  scale_x_continuous(expand = c(-1, 107))+
  geom_vline(xintercept=42,linetype="dashed")+
  scale_fill_cyclical(
    breaks = c("tBRCAm (n=863)","Non-BRCA HRRm (n=107)","Non-HRRm (n=868)"),
    #labels = c(PAOLA = "PAOLA", SCR = "Screening"),
    values = c("dodgerblue2","orange","forestgreen"),
    name = "Status",
    guide = "legend"
  ) +
  theme_ridges(grid = TRUE)+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  guides(fill=guide_legend(title=""))+
  theme(strip.text.x = element_text(size=20, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=2, linetype="solid"))+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = 'top')





ggplot(all_HRD_plots,aes(y=Status))+
  geom_density_ridges(scale= 4,
                      aes(x = HRD_Score, fill =Status),
                      alpha = .8, color = "black", from = 0, to=100) +
  labs(
    x = "HRD Score",
    y = "",
    title = "",
    subtitle = "",
    caption = ""
  )+
  facet_wrap(~Study,scales="free_x",ncol=1) +
  scale_y_discrete(expand = expand_scale(mult = c(0.01, 1.8)))+ 
  #scale_y_discrete(expand = c(0.1, -0.1))+
  scale_x_continuous(expand = c(-1, 104))+
  scale_fill_cyclical(
    breaks = c("tBRCAm","HRRm and non-BRCAm","Non-HRRm"),
    #labels = c(PAOLA = "PAOLA", SCR = "Screening"),
    values = c("dodgerblue2","orange","forestgreen"),
    name = "Status",
    guide = "legend"
  ) +
  theme_ridges(grid = TRUE)+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  guides(fill=guide_legend(title=""))+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=2, linetype="solid"))+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')

all_HRD_plots_BRCA<-subset(all_HRD_plots,all_HRD_plots$Status_3 =="tBRCAm")
all_HRD_plots_BRCA_LOH<-all_HRD_plots_BRCA[!all_HRD_plots_BRCA$LOH =="Unknown",]
all_HRD_plots_BRCA_LOH<-all_HRD_plots_BRCA_LOH[!is.na(all_HRD_plots_BRCA_LOH$LOH),]

t.test(all_HRD_plots_BRCA_LOH$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Biallelic Inactivation"],all_HRD_plots$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Heterozygous"])

t.test(all_HRD_plots_BRCA_LOH$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Biallelic Inactivation" & all_HRD_plots_BRCA_LOH$Status_2 =="BRCA1m"],all_HRD_plots$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Heterozygous" & all_HRD_plots_BRCA_LOH$Status_2 =="BRCA1m"])
t.test(all_HRD_plots_BRCA_LOH$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Biallelic Inactivation" & all_HRD_plots_BRCA_LOH$Status_2 =="BRCA2m"],all_HRD_plots$HRD_Score[all_HRD_plots_BRCA_LOH$LOH =="Heterozygous" & all_HRD_plots_BRCA_LOH$Status_2 =="BRCA2m"])


ggplot(all_HRD_plots_BRCA_LOH, aes(x = LOH, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=LOH),outlier.shape = NA)+
  xlab("LOH status")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  facet_grid(.~ Status_2,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  scale_fill_manual(values=c("Biallelic Inactivation"="dodgerblue","Heterozygous" ="grey49"))+
    theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("GIS")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

all_HRD_plots_BRCA_LOH1<-all_HRD_plots_BRCA_LOH[!all_HRD_plots_BRCA_LOH$S_G_Status =="s/g BRCA not determined",]
ggplot(all_HRD_plots_BRCA_LOH1, aes(x = LOH, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=LOH),outlier.shape = NA)+
  xlab("LOH status")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  facet_grid(.~ S_G_Status,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  scale_fill_manual(values=c("Biallelic Inactivation"="dodgerblue","Heterozygous" ="grey49"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("GIS")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


SOLOs<-c("SOLO1","SOLO2")  
all_HRD_nsolo<-all_HRD[!all_HRD$Study %in% SOLOs,]

  
take_out<-c("BRCA1m ; BRCA2m","Test Failed/Cancelled","BRCA1m;tBRCA2m")  
all_HRD_plots_1<-all_HRD_plots[!all_HRD_plots$Status_2 %in% take_out,]

  
  ####boxplots
  #####all
####tBRCAm vs nontBRCAm
ggplot(all_HRD_plots, aes(x = Status_3, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=Status_3),outlier.shape = NA)+
  xlab("Subgroup")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  scale_fill_manual(values=c("tBRCAm"="dodgerblue2","Non-tBRCAm" ="limegreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

t.test(all_HRD_plots$HRD_Score[all_HRD_plots$Status_3 =="tBRCAm"],all_HRD_plots$HRD_Score[all_HRD_plots$Status_3 =="Non-tBRCAm"])







ggplot(all_HRD_plots, aes(x = Status_3, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=Status_3),outlier.shape = NA)+
  xlab("Subgroup")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  scale_fill_manual(values=c("Biallelic Inactivation"="dodgerblue2","Heterozygous" ="grey49"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))



  ggplot(all_HRD_plots_1, aes(x = Status_2, y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=Status_2),outlier.shape = NA)+
    xlab("Subgroup")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Study,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+
    theme(legend.title=element_text(size=14,face="bold"))+
    theme(legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=20,face="bold"))+
    theme(legend.text=element_text(size=20))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))

  
  all_HRD_plots_1$Status_4<-case_when(all_HRD_plots_1$Status_2 =="BRCA1m"~"BRCA1m",
                                      all_HRD_plots_1$Status_2 =="BRCA2m"~"BRCA2m",
                                      all_HRD_plots_1$Status_3 =="Non-tBRCAm"~"Non-tBRCAm")
  
  ggplot(all_HRD_plots_1, aes(x = Status_4, y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=Status_4),outlier.shape = NA)+
    xlab("Subgroup")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Study,scales="free")+
    facet_grid(.~ Population,scales="free")+
    
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2","Non-tBRCAm" ="limegreen"))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"))+
    theme(legend.title=element_text(size=14,face="bold"))+
    theme(legend.text=element_text(size=14))+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=20,face="bold"))+
    theme(legend.text=element_text(size=20))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))
  
  
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="Non-tBRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="BRCA1m"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="Non-tBRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="BRCA2m"])  
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="BRCA1m"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_4 =="BRCA2m"])  
  
  
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status =="HRRm and non-BRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status =="tBRCAm"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status =="HRRm and non-BRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status =="Non-HRRm"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="HRRm and non-BRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="HRRm and non-BRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="HRRm and non-BRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="Non-HRRm"])
  
  
    
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Population =="CR/PR after 1L Platinum Based Chemotherapy"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Population =="CR/PR after 1L Platinum Based Chemotherapy"])
  
  
    ##### 6 studies
    ggplot(all_HRD_plots_1, aes(x = Status_2, y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=Status_2),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    facet_grid(.~ Study,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=20, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=8),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=20),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=20))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))

  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="PAOLA"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="PAOLA"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="OPINION"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="OPINION"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="LIGHT"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="LIGHT"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="Study 19"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="Study 19"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="SOLO1"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="SOLO1"])
  t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Study =="SOLO2"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Study =="SOLO2"])
  

  ##### Population
  ggplot(all_HRD_plots_1, aes(x = Status_2, y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=Status_2),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    facet_grid(.~ Population,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=23, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=20),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=20),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=20))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))  

t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Population =="CR/PR after 1L Platinum Based Chemotherapy"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA1m" & all_HRD_plots_1$Population =="PSR"])
t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Population =="CR/PR after 1L Platinum Based Chemotherapy"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="BRCA2m" & all_HRD_plots_1$Population =="PSR"])
t.test(all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="HRRm and non-BRCAm" & all_HRD_plots_1$Population =="CR/PR after 1L Platinum Based Chemotherapy"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$Status_2 =="HRRm and non-BRCAm" & all_HRD_plots_1$Population =="PSR"])

    
  
###G/S BRCAm
  all_HRD_plots_2<-all_HRD_plots_1[!is.na(all_HRD_plots_1$S_G_Status),]  
  all_HRD_plots_2$S_G_Status <- factor(all_HRD_plots_2$S_G_Status, levels=c( "gBRCAm","sBRCAm","s/g BRCA not determined","HRRm and non-BRCAm","Non-HRRm")) #reordering zygosity variables to appear in this order in the barplot
  
    
  ggplot(all_HRD_plots_2, aes(x = S_G_Status , y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=S_G_Status),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Population,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("gBRCAm"="blueviolet","sBRCAm"="cyan3",`s/g BRCA not determined`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=20),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=20),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=20,face="bold"))+
    theme(legend.text=element_text(size=20))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))

  t.test(all_HRD_plots_2$HRD_Score[all_HRD_plots_2$S_G_Status =="gBRCAm"],all_HRD_plots_1$HRD_Score[all_HRD_plots_1$S_G_Status =="sBRCAm"])
  
  
  ###G/S BRCAm
  all_HRD_plots_3<-all_HRD_plots_2[!is.na(all_HRD_plots_2$S_G_Status_BRCA1_2),]  
  all_HRD_plots_3$S_G_Status_BRCA1_2 <- factor(all_HRD_plots_3$S_G_Status_BRCA1_2, levels=c( "gBRCA1m","gBRCA2m","sBRCA1m","sBRCA2m","s/g BRCA not determined BRCA1m","s/g BRCA not determined BRCA2m","HRRm and non-BRCAm","Non-HRRm")) #reordering zygosity variables to appear in this order in the barplot
  
  all_HRD_plots_3$S_G_Status_BRCA1_3<-case_when(all_HRD_plots_3$S_G_Status_BRCA1_2 =="gBRCA1m" ~ "gBRCA1m",
                                                all_HRD_plots_3$S_G_Status_BRCA1_2 =="gBRCA2m" ~ "gBRCA2m",
                                                all_HRD_plots_3$S_G_Status_BRCA1_2 =="sBRCA1m" ~ "sBRCA1m",
                                                all_HRD_plots_3$S_G_Status_BRCA1_2 =="sBRCA2m" ~ "sBRCA2m",
                                                all_HRD_plots_3$S_G_Status_BRCA1_2 =="s/g BRCA not determined BRCA1m" ~ "s/g BRCA not determined BRCA1m",
                                                all_HRD_plots_3$S_G_Status_BRCA1_2 =="s/g BRCA not determined BRCA2m" ~ "s/g BRCA not determined BRCA2m",
                                                all_HRD_plots_3$Status_4 =="Non-tBRCAm" ~ "Non-tBRCAm")
  
  
  
  ggplot(all_HRD_plots_3, aes(x = S_G_Status_BRCA1_2 , y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=S_G_Status_BRCA1_2),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Study,scales="free")
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=12, face= "bold"),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=16),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=16))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))
  
  all_HRD_plots_3$S_G_Status_BRCA1_3 <- factor(all_HRD_plots_3$S_G_Status_BRCA1_3, levels=c( "gBRCA1m","gBRCA2m","sBRCA1m","sBRCA2m","s/g BRCA not determined BRCA1m","s/g BRCA not determined BRCA2m","Non-tBRCAm")) #reordering zygosity variables to appear in this order in the barplot
  
  
  ggplot(all_HRD_plots_3, aes(x = S_G_Status_BRCA1_3 , y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=S_G_Status_BRCA1_3),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Study,scales="free")
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3","s/g BRCA not determined BRCA1m"= "grey50","s/g BRCA not determined BRCA2m"= "grey50","Non-tBRCAm" ="limegreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=12, face= "bold"),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=16),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=16))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))
  
  
  
all_HRD_plots_3_BRCA<-all_HRD_plots_3[all_HRD_plots_3$Status =="tBRCAm",]
all_HRD_plots_3_BRCA<-all_HRD_plots_3_BRCA[!all_HRD_plots_3_BRCA$S_G_Status =="s/g BRCA not determined",]
all_HRD_plots_3_BRCA<-all_HRD_plots_3_BRCA[!all_HRD_plots_3_BRCA$LOH =="Unknown",]

  ggplot(all_HRD_plots_3_BRCA, aes(x = paste(S_G_Status_BRCA1_2,LOH) , y = HRD_Score),color="black")+
    geom_boxplot(aes(fill=LOH),outlier.shape = NA)+
    xlab("")+ 
    #coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    facet_wrap(.~ S_G_Status_BRCA1_2,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
    theme(panel.spacing = unit(6, "lines"))+
    #theme_bw()+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=12, face= "bold"),axis.title.x=element_text(size=12))+
    theme(axis.text.y=element_text(size=16),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=16))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))
  
  
  ###violin plot subtypes of variants
  all_HRD_plots_4<-all_HRD_plots[!is.na(all_HRD_plots$Variant_Subtype),]
  all_HRD_plots_4_BRCAm<-subset(all_HRD_plots_4,all_HRD_plots_4$Status =="tBRCAm")
  
  all_HRD_plots_4_BRCAm$Variant_Subtype <- factor(all_HRD_plots_4_BRCAm$Variant_Subtype, levels=c( "Frameshift/Indel","Nonsense","Loss/Rearrangement","Splice","Missense")) #reordering zygosity variables to appear in this order in the barplot
  #all_HRD_plots_4_BRCAm<-all_HRD_plots_4[!all_HRD_plots_4$Status_2 %in% take_out,]
  
  ggplot(all_HRD_plots_4_BRCAm, aes(x = Variant_Subtype , y = HRD_Score),color="black")+
    geom_violin(aes(fill=Variant_Subtype),size=1)+
    xlab("")+ 
    coord_flip()+
    geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
    geom_hline(yintercept=42,linetype="dashed")+
    #facet_grid(.~ Status_3,scales="free")+
    scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
    #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
    theme_bw()+
    theme(strip.text.x = element_text(size=15, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=1.5, linetype="solid"))+
    theme(axis.text.x=element_text(size=14),axis.title.x=element_text(size=20,face="bold"))+
    theme(axis.text.y=element_text(size=16),axis.title.y =element_text(size=20,face="bold"))+
    theme(legend.title=element_text(size=19,face="bold"))+
    theme(legend.text=element_text(size=16))+
    ylab("HRD Score")+
    theme(legend.position = 'top')+
    guides(fill=guide_legend(title=""))
  


  ###LOH data

  all_HRD_zyg<-all_HRD[!all_HRD$LOH =="Unknown",]
  all_HRD_zyg_1<-all_HRD_zyg[!all_HRD_zyg$LOH =="",]
  all_HRD_zyg_2<-all_HRD_zyg_1[!is.na(all_HRD_zyg_1$LOH),]
  all_HRD_zyg_2$GENE_2<-ifelse(all_HRD_zyg_2$GENE =="BRCA1;BRCA1", "BRCA1",
                               ifelse(all_HRD_zyg_2$GENE =="BRCA2;BRCA2", "BRCA2",
                                      ifelse(all_HRD_zyg_2$GENE =="CDK12;CDK12", "CDK12",as.character(all_HRD_zyg_2$GENE))))
  
  
  
  all_1<-ddply(all_HRD_zyg_2, .(all_HRD_zyg_2$GENE_2,all_HRD_zyg_2$LOH),nrow)
  colnames(all_1)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes
  
  dd4=as.data.table(all_1) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
    dcast(Gene~Zygosity,value.var="Count")
  

  ### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
  dd4[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
  dd4[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases

  
  
  dd_df4<-data.frame(dd4) ##reformatting the data table (dd3) into a dataframe
  dd_df4_melt<-melt(dd_df4, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
  dd_df4_melt_ord<-dd_df4_melt[order(dd_df4_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
  dd_df4_melt_ord_frac<-subset(dd_df4_melt_ord,dd_df4_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
  dd_df4_melt_ord_count<-dd_df4_melt_ord[!dd_df4_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
  
  dd_df4_melt_ord_frac$variable<-gsub("_frac","",dd_df4_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
  dd_df4_melt_ord_frac$key<-paste(dd_df4_melt_ord_frac$Gene,dd_df4_melt_ord_frac$variable,sep="_")
  dd_df4_melt_ord_count$key<-paste(dd_df4_melt_ord_count$Gene,dd_df4_melt_ord_count$variable,sep="_")
  dd_df4_melt_ord_frac_merged<-merge(dd_df4_melt_ord_frac,dd_df4_melt_ord_count,by="key")
  dd_df4_melt_ord_frac<-dd_df4_melt_ord_frac_merged[,c(2:4,7)]
  
  dd_df4_melt_ord_frac_new<-dd_df4_melt_ord_frac[!is.na(dd_df4_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting
  
  #dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting
  
  
  colnames(dd_df4_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
  dd_df4_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df4_melt_ord_frac_new$Zygosity)
  
  dd_df4_melt_ord_frac_new$Zygosity <- factor(dd_df4_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
  dd_df4_melt_ord_frac_new$Gene <- factor(dd_df4_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot
  
  dd_df4_melt_ord_frac_new$Status<-ifelse(dd_df4_melt_ord_frac_new$Gene %in% BRCA,"tBRCAm","HRRm and non-BRCAm")
  dd_df4_melt_ord_frac_new$Status <- gsub("HRRm and non-BRCAm","Non-BRCA HRRm", dd_df4_melt_ord_frac_new$Status)
  
  dd_df4_melt_ord_frac_new$Status <- factor(dd_df4_melt_ord_frac_new$Status, levels=c("tBRCAm","Non-BRCA HRRm")) #reordering zygosity variables to appear in this order in the barplot
  
  #geom_text(aes(label = paste0(Freq,"\n (", Perc,"%",")")), size = 9,color="white", position = position_stack(vjust = 0.45))+
    
  
  ggplot(dd_df4_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity))+ 
    geom_bar(stat = "identity",color="black",size=0.1,width =0.7)+
    xlab("\nGene")+ 
    ylab("Percent (%)\n")+
    #geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))
    scale_fill_manual(values=c("Biallelic Inactivation"="dodgerblue", "Heterozygous"="maroon4"), breaks=c("Biallelic Inactivation","Heterozygous"))+
    geom_text(aes(label = paste0(Count,"\n (", Value,"%",")")), size = 6.9,color="white", position = position_stack(vjust = 0.45))+
    facet_grid(.~ Status,scales="free",space="free")+
    theme_bw()+
    theme(panel.spacing = unit(3, "lines"))+
    guides(color=guide_legend(title="Zygosity"))+
    theme(strip.text.x = element_text(size=25, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=3, linetype="solid"))+
    theme(axis.text=element_text(size=24),axis.title.y=element_text(size=24,face="bold"),axis.title.x=element_blank())+
    theme(legend.title=element_text(size=24,face="bold"))+
    theme(legend.text=element_text(size=24))+
    theme(legend.position = 'top')
 
    
    
 #####plotting count rather than proportion
  
  ggplot(dd_df4_melt_ord_frac_new, aes(x = Gene, y = Count, fill = Zygosity)) + 
    geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
    xlab("\nGene")+ 
    ylab("Number of HRRm cases with evaluable zygosity\n")+
    #geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))
    scale_fill_manual(values=c("Biallelic Inactivation"="dodgerblue", "Heterozygous"="maroon4"), breaks=c("Biallelic Inactivation","Heterozygous"))+
    facet_grid(.~ Status,scales="free",space="free")+
    theme_bw()+
    theme(panel.spacing = unit(3, "lines"))+
    guides(color=guide_legend(title="Zygosity"))+
    theme(strip.text.x = element_text(size=25, color="Black",
                                      face="bold.italic"))+
    theme(strip.background = element_rect(colour="black", fill="white", 
                                          size=3, linetype="solid"))+
    theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
    theme(legend.title=element_text(size=18,face="bold"))+
    theme(legend.text=element_text(size=18))+
    theme(legend.position = 'top')
  
  write.csv(dd4,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/gsLOH_no_unknowns_1plot.txt",row.names=FALSE) 

  
#####now do study by study
all_HRD_zyg_2_PAOLA<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="PAOLA")
all_HRD_zyg_2_OP<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="OPINION")
all_HRD_zyg_2_L<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="LIGHT")
all_HRD_zyg_2_19<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="Study 19")
all_HRD_zyg_2_S1<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="SOLO1")
all_HRD_zyg_2_S2<-subset(all_HRD_zyg_2,all_HRD_zyg_2$Study =="SOLO2")


###PAOLA
all_1_P<-ddply(all_HRD_zyg_2_PAOLA, .(all_HRD_zyg_2_PAOLA$GENE_2,all_HRD_zyg_2_PAOLA$LOH),nrow)
colnames(all_1_P)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd5=as.data.table(all_1_P) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd5[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd5[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df5<-data.frame(dd5) ##reformatting the data table (dd3) into a dataframe
dd_df5_melt<-melt(dd_df5, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df5_melt_ord<-dd_df5_melt[order(dd_df5_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df5_melt_ord_frac<-subset(dd_df5_melt_ord,dd_df5_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df5_melt_ord_count<-dd_df5_melt_ord[!dd_df5_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df5_melt_ord_frac$variable<-gsub("_frac","",dd_df5_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df5_melt_ord_frac$key<-paste(dd_df5_melt_ord_frac$Gene,dd_df5_melt_ord_frac$variable,sep="_")
dd_df5_melt_ord_count$key<-paste(dd_df5_melt_ord_count$Gene,dd_df5_melt_ord_count$variable,sep="_")
dd_df5_melt_ord_frac_merged<-merge(dd_df5_melt_ord_frac,dd_df5_melt_ord_count,by="key")
dd_df5_melt_ord_frac<-dd_df5_melt_ord_frac_merged[,c(2:4,7)]

dd_df5_melt_ord_frac_new<-dd_df5_melt_ord_frac[!is.na(dd_df5_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df5_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df5_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df5_melt_ord_frac_new$Zygosity)

dd_df5_melt_ord_frac_new$Zygosity <- factor(dd_df5_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df5_melt_ord_frac_new$Gene <- factor(dd_df5_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","RAD51C","RAD51D","BRIP1","ATM","PALB2","CHEK2","RAD51B","FANCL","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df5_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title.y=element_text(size=18,face="bold"),axis.title.x=element_blank())+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')

axis.title.x=element_blank()
###OPINION
all_1_O<-ddply(all_HRD_zyg_2_OP, .(all_HRD_zyg_2_OP$GENE_2,all_HRD_zyg_2_OP$LOH),nrow)
colnames(all_1_O)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd6=as.data.table(all_1_O) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd6[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd6[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df6<-data.frame(dd6) ##reformatting the data table (dd3) into a dataframe
dd_df6_melt<-melt(dd_df6, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df6_melt_ord<-dd_df6_melt[order(dd_df6_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df6_melt_ord_frac<-subset(dd_df6_melt_ord,dd_df6_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df6_melt_ord_count<-dd_df6_melt_ord[!dd_df6_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df6_melt_ord_frac$variable<-gsub("_frac","",dd_df6_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df6_melt_ord_frac$key<-paste(dd_df6_melt_ord_frac$Gene,dd_df6_melt_ord_frac$variable,sep="_")
dd_df6_melt_ord_count$key<-paste(dd_df6_melt_ord_count$Gene,dd_df6_melt_ord_count$variable,sep="_")
dd_df6_melt_ord_frac_merged<-merge(dd_df6_melt_ord_frac,dd_df6_melt_ord_count,by="key")
dd_df6_melt_ord_frac<-dd_df6_melt_ord_frac_merged[,c(2:4,7)]

dd_df6_melt_ord_frac_new<-dd_df6_melt_ord_frac[!is.na(dd_df6_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df6_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df6_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df6_melt_ord_frac_new$Zygosity)

dd_df6_melt_ord_frac_new$Zygosity <- factor(dd_df6_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df6_melt_ord_frac_new$Gene <- factor(dd_df6_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df6_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')
 


###LIGHT
all_1_L<-ddply(all_HRD_zyg_2_L, .(all_HRD_zyg_2_L$GENE_2,all_HRD_zyg_2_L$LOH),nrow)
colnames(all_1_L)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd7=as.data.table(all_1_L) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd7[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd7[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df7<-data.frame(dd7) ##reformatting the data table (dd3) into a dataframe
dd_df7_melt<-melt(dd_df7, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df7_melt_ord<-dd_df7_melt[order(dd_df7_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df7_melt_ord_frac<-subset(dd_df7_melt_ord,dd_df7_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df7_melt_ord_count<-dd_df7_melt_ord[!dd_df7_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df7_melt_ord_frac$variable<-gsub("_frac","",dd_df7_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df7_melt_ord_frac$key<-paste(dd_df7_melt_ord_frac$Gene,dd_df7_melt_ord_frac$variable,sep="_")
dd_df7_melt_ord_count$key<-paste(dd_df7_melt_ord_count$Gene,dd_df7_melt_ord_count$variable,sep="_")
dd_df7_melt_ord_frac_merged<-merge(dd_df7_melt_ord_frac,dd_df7_melt_ord_count,by="key")
dd_df7_melt_ord_frac<-dd_df7_melt_ord_frac_merged[,c(2:4,7)]

dd_df7_melt_ord_frac_new<-dd_df7_melt_ord_frac[!is.na(dd_df7_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df7_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df7_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df7_melt_ord_frac_new$Zygosity)

dd_df7_melt_ord_frac_new$Zygosity <- factor(dd_df7_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df7_melt_ord_frac_new$Gene <- factor(dd_df7_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df7_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')


####Study 19
all_1_19<-ddply(all_HRD_zyg_2_19, .(all_HRD_zyg_2_19$GENE_2,all_HRD_zyg_2_19$LOH),nrow)
colnames(all_1_19)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd8=as.data.table(all_1_19) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd8[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd8[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df8<-data.frame(dd8) ##reformatting the data table (dd3) into a dataframe
dd_df8_melt<-melt(dd_df8, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df8_melt_ord<-dd_df8_melt[order(dd_df8_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df8_melt_ord_frac<-subset(dd_df8_melt_ord,dd_df8_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df8_melt_ord_count<-dd_df8_melt_ord[!dd_df8_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df8_melt_ord_frac$variable<-gsub("_frac","",dd_df8_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df8_melt_ord_frac$key<-paste(dd_df8_melt_ord_frac$Gene,dd_df8_melt_ord_frac$variable,sep="_")
dd_df8_melt_ord_count$key<-paste(dd_df8_melt_ord_count$Gene,dd_df8_melt_ord_count$variable,sep="_")
dd_df8_melt_ord_frac_merged<-merge(dd_df8_melt_ord_frac,dd_df8_melt_ord_count,by="key")
dd_df8_melt_ord_frac<-dd_df8_melt_ord_frac_merged[,c(2:4,7)]

dd_df8_melt_ord_frac_new<-dd_df8_melt_ord_frac[!is.na(dd_df8_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df8_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df8_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df8_melt_ord_frac_new$Zygosity)

dd_df8_melt_ord_frac_new$Zygosity <- factor(dd_df8_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df8_melt_ord_frac_new$Gene <- factor(dd_df8_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df8_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')


####SOLO1
all_1_S1<-ddply(all_HRD_zyg_2_S1, .(all_HRD_zyg_2_S1$GENE_2,all_HRD_zyg_2_S1$LOH),nrow)
colnames(all_1_S1)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd9=as.data.table(all_1_S1) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd9[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd9[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df9<-data.frame(dd9) ##reformatting the data table (dd3) into a dataframe
dd_df9_melt<-melt(dd_df9, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df9_melt_ord<-dd_df9_melt[order(dd_df9_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df9_melt_ord_frac<-subset(dd_df9_melt_ord,dd_df9_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df9_melt_ord_count<-dd_df9_melt_ord[!dd_df9_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df9_melt_ord_frac$variable<-gsub("_frac","",dd_df9_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df9_melt_ord_frac$key<-paste(dd_df9_melt_ord_frac$Gene,dd_df9_melt_ord_frac$variable,sep="_")
dd_df9_melt_ord_count$key<-paste(dd_df9_melt_ord_count$Gene,dd_df9_melt_ord_count$variable,sep="_")
dd_df9_melt_ord_frac_merged<-merge(dd_df9_melt_ord_frac,dd_df9_melt_ord_count,by="key")
dd_df9_melt_ord_frac<-dd_df9_melt_ord_frac_merged[,c(2:4,7)]

dd_df9_melt_ord_frac_new<-dd_df9_melt_ord_frac[!is.na(dd_df9_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df9_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df9_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df9_melt_ord_frac_new$Zygosity)

dd_df9_melt_ord_frac_new$Zygosity <- factor(dd_df9_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df9_melt_ord_frac_new$Gene <- factor(dd_df9_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df9_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')

###SOLO2   

all_1_S2<-ddply(all_HRD_zyg_2_S2, .(all_HRD_zyg_2_S2$GENE_2,all_HRD_zyg_2_S2$LOH),nrow)
colnames(all_1_S2)<-c("Gene","Zygosity","Count") ##renaming column names in object HRR_genes_Z_1 for formatting purposes

dd10=as.data.table(all_1_S2) %>%  # creatinga  data table (dd3) that splits object (HRR_genes_Z_1) into separate columns specific to each of the 4 zygosity variables
  dcast(Gene~Zygosity,value.var="Count")


### adding additional columns to data table (dd3) where the proportion of Biallelic Inactivation,Suspected Biallelic Inactivation,Heterozygous and unknown variables are calculated based on the total count of each of these cases
dd10[,Biallelic.Inactivation_frac:=round((`Biallelic Inactivation`/sum(na.rm =T, `Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ### (Biallelic_Inactivation_frac) is the fraction of Biallelic Inactivation cases relative to the sum of all zygosity cases
dd10[,Heterozygous_frac:=round((Heterozygous/sum(na.rm =T,`Biallelic Inactivation`,Heterozygous))*100,digits=2),Gene] ##(Heterozygous_frac) is the fraction of Heterozygous cases relative to the sum of all zygosity cases



dd_df10<-data.frame(dd10) ##reformatting the data table (dd3) into a dataframe
dd_df10_melt<-melt(dd_df10, id="Gene") # reformatting the (dd_df3) object into a molten dataframe. This is for ease of plotting purposes  
dd_df10_melt_ord<-dd_df10_melt[order(dd_df10_melt$Gene),] ##Ordering the molten dataframe (dd_df3_melt) by Gene,expressed as new object (dd_df3_melt_ord)
dd_df10_melt_ord_frac<-subset(dd_df10_melt_ord,dd_df10_melt_ord$variable %like% "frac") ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df10_melt_ord_count<-dd_df10_melt_ord[!dd_df10_melt_ord$variable %like% "frac",] ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot

dd_df10_melt_ord_frac$variable<-gsub("_frac","",dd_df10_melt_ord_frac$variable) ###filtering the object (dd_df3_melt_ord) for just the fractions for each zygosity variable as this is what we will plot
dd_df10_melt_ord_frac$key<-paste(dd_df10_melt_ord_frac$Gene,dd_df10_melt_ord_frac$variable,sep="_")
dd_df10_melt_ord_count$key<-paste(dd_df10_melt_ord_count$Gene,dd_df10_melt_ord_count$variable,sep="_")
dd_df10_melt_ord_frac_merged<-merge(dd_df10_melt_ord_frac,dd_df10_melt_ord_count,by="key")
dd_df10_melt_ord_frac<-dd_df10_melt_ord_frac_merged[,c(2:4,7)]

dd_df10_melt_ord_frac_new<-dd_df10_melt_ord_frac[!is.na(dd_df10_melt_ord_frac$value.x),] #creating a new object (dd_df3_melt_ord_frac_new) where missing values (NA) in the "value" column are removed from (dd_df3_melt_ord_frac) to help with plotting

#dd_df3_melt_ord_frac_new$variable<-gsub("_frac","",dd_df3_melt_ord_frac_new$variable) ##removing the "frac" from the name of each variable in object (dd_df3_melt_ord_frac_new) to help with plotting


colnames(dd_df10_melt_ord_frac_new)<-c("Gene","Zygosity","Value","Count") #renaming columns to help with plotting
dd_df10_melt_ord_frac_new$Zygosity<-gsub("Biallelic.Inactivation","Biallelic Inactivation",dd_df10_melt_ord_frac_new$Zygosity)

dd_df10_melt_ord_frac_new$Zygosity <- factor(dd_df10_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
dd_df10_melt_ord_frac_new$Gene <- factor(dd_df10_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot


ggplot(dd_df10_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_classic()+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=25, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=10, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')



#####
dd_df5_melt_ord_frac_new$Study <-"PAOLA"
dd_df6_melt_ord_frac_new$Study <-"OPINION"
dd_df7_melt_ord_frac_new$Study<-"LIGHT"
dd_df8_melt_ord_frac_new$Study <-"Study 19"
dd_df9_melt_ord_frac_new$Study <-"SOLO1"
dd_df10_melt_ord_frac_new$Study <-"SOLO2"


dd_df_melt_ord_frac_new<-rbind(dd_df5_melt_ord_frac_new,dd_df6_melt_ord_frac_new,dd_df7_melt_ord_frac_new,dd_df8_melt_ord_frac_new,dd_df9_melt_ord_frac_new,dd_df10_melt_ord_frac_new)

write.csv(dd_df_melt_ord_frac_new,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/dd_df_melt_ord_frac_new_newColours.csv",row.names=FALSE)

dd_df_melt_ord_frac_new<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/dd_df_melt_ord_frac_new_newColours.csv")

dd_df_melt_ord_frac_new$Zygosity <- factor(dd_df_melt_ord_frac_new$Zygosity, levels=c("Heterozygous","Biallelic Inactivation")) #reordering zygosity variables to appear in this order in the barplot
#dd_df_melt_ord_frac_new$Gene <- factor(dd_df_melt_ord_frac_new$Gene, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot
dd_df_melt_ord_frac_new$Gene <- factor(dd_df_melt_ord_frac_new$Gene, levels=c("PPP2R2A","BARD1","RAD54L","CHEK2","FANCL","RAD51B","PALB2","ATM","RAD51D","RAD51C","BRIP1","CDK12","BRCA2","BRCA1")) #reordering HRR Gene variables to appear in this order in the barplot

dd_df_melt_ord_frac_new$Study <- factor(dd_df_melt_ord_frac_new$Study, levels=c("PAOLA","OPINION","LIGHT","Study 19","SOLO1","SOLO2")) #reordering HRR Gene variables to appear in this order in the barplot

dd_df_melt_ord_frac_new$Status<-ifelse(dd_df_melt_ord_frac_new$Gene %in% BRCA,"tBRCAm","HRRm and non-BRCAm")
dd_df_melt_ord_frac_new$Status <- factor(dd_df_melt_ord_frac_new$Status, levels=c("tBRCAm","HRRm and non-BRCAm")) #reordering zygosity variables to appear in this order in the barplot


ggplot(dd_df_melt_ord_frac_new, aes(x = Gene, y = Value, fill = Zygosity)) + 
  geom_bar(stat = "identity",color="black",size=0.1,width =0.7,alpha=0.8)+
  xlab("\nGene")+ 
  ylab("Percent (%)\n")+
  geom_text(aes(label = Count), size = 5,color="white", position = position_stack(vjust = 0.5))+
  #geom_text(aes(label = paste0(Count," (", Value,"%",")")), size = 4.9,color="white", position = position_stack(vjust = 0.45))+
  coord_flip()+
  #facet_grid(Status ~ Study,scales="free", space="free")+
  facet_grid(.~ Study,scales="free")+
  scale_fill_manual(values=c("Biallelic Inactivation"="blue4", "Heterozygous"="grey"), breaks=c("Biallelic Inactivation","Heterozygous"))+
  theme_bw()+
  theme(panel.spacing = unit(2, "lines"))+
  theme(legend.position = "top")+
  theme(legend.text=element_text(size=14))+
  guides(color=guide_legend(title="Zygosity"))+
  theme(strip.text.x = element_text(size=22, color="Black",
                                    face="bold.italic"))+
  #theme(strip.text.y = element_text(size=22, color="Black",
   #                                 face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=3, linetype="solid"))+
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18,face="bold"))+
  theme(legend.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=18))+
  theme(legend.position = 'top')



###HRRm
all_HRD_HRR<-subset(all_HRD, all_HRD$Status =="HRRm and non-BRCAm")
all_HRD_HRR<-subset(all_HRD, all_HRD$Status =="HRRm and non-BRCAm" | all_HRD$Status =="tBRCAm")

all_HRD_plots_HRR<-all_HRD_plots[all_HRD_plots$Status =="Non-BRCA HRRm",]
all_HRD_plots_HRR<-all_HRD_plots[all_HRD_plots$Status =="Non-BRCA HRRm" | all_HRD_plots$Status =="tBRCAm",]

ggplot(all_HRD_plots_HRR, aes(x = Study, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=Study),outlier.shape = NA)+
  xlab("Study")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(legend.position = "top")+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="Study 19"])


t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="Study 19"])

t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="PAOLA"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="OPINION"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"])
t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="LIGHT"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Study =="Study 19"])


t.test(all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Population =="CR/PR after 1L Platinum Based Chemotherapy"],all_HRD_plots_HRR$HRD_Score[all_HRD_plots_HRR$Population =="PSR"])


tiff("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Figures/Plot10_density_all_studiesHRRm.tiff", units="in", width=5, height=5, res=200)
ggplot(all_HRD_plots_HRR, aes(x=HRD_Score))+
  geom_density(aes(fill=Study),color="black",alpha=0.3)+
  xlab("HRD Score")+
  ylab("Normalised Density")+
  theme_bw()+
  scale_x_continuous(limits = c(0, 101), breaks = seq(0, 101, by = 25))+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=9,face="bold"))+
  theme(legend.title=element_text(size=9,face="bold"))+
  theme(legend.text=element_text(size=9))+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))
dev.off()



ggplot(all_HRD_plots_HRR, aes(x = Population, y = HRD_Score),color="black")+
  geom_boxplot(aes(fill=Population),outlier.shape = NA)+
  xlab("")+ 
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(legend.position = "top")+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


tiff("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Figures/Plot12_density_all_studiesHRRm_1L_PSR.tiff", units="in", width=5, height=5, res=200)
ggplot(all_HRD_plots_HRR, aes(x=HRD_Score))+
  geom_density(aes(fill=Population),color="black",alpha=0.3)+
  xlab("HRD Score")+
  ylab("Normalised Density")+
  theme_bw()+
  scale_x_continuous(limits = c(0, 101), breaks = seq(0, 101, by = 25))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=9,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))
dev.off()




all_HRD_plots_HRR$GENE_2<-ifelse(all_HRD_plots_HRR$GENE =="CDK12;CDK12","CDK12",
                                 ifelse(all_HRD_plots_HRR$GENE %like% ";", "Co-occurring",as.character(all_HRD_plots_HRR$GENE)))


all_HRD_plots_HRR_1<-all_HRD_plots_HRR[!all_HRD_plots_HRR$GENE_2 =="Co-occurring",]

all_HRD_plots_HRR_1$GENE_2 <- factor(all_HRD_plots_HRR_1$GENE_2, levels=c("BRCA1","BRCA2","CDK12","BRIP1","RAD51C","RAD51D","ATM","PALB2","RAD51B","FANCL","CHEK2","RAD54L","BARD1","PPP2R2A")) #reordering HRR Gene variables to appear in this order in the barplot
names(all_HRD_plots_HRR_1)[names(all_HRD_plots_HRR_1) == 'LOH'] <- 'Zygosity'



ggplot(all_HRD_plots_HRR_1, aes(x = GENE_2, y = HRD_Score),color="black")+
  geom_boxplot(outlier.shape = NA)+
  xlab("GENE")+ 
  #coord_flip()+
  geom_jitter(aes(color= Zygosity),shape=16, position=position_jitter(0.2),size=3)+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study.x,scales="free")+
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


ggplot(all_HRD_plots_HRR_1, aes(x = GENE_2, y = HRD_Score),color="black")+
  geom_boxplot(outlier.shape = NA)+
  xlab("GENE")+ 
  coord_flip()+
  geom_jitter(aes(color= Zygosity),shape=16, position=position_jitter(0.2),size=3)+
  geom_hline(yintercept=42,linetype="dashed")+
  facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=22, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=22,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


###

sub<-c("BRIP1","RAD51C","RAD51D","PALB2")

all_HRD_plots_HRR_1_sub<-subset(all_HRD_plots_HRR_1,all_HRD_plots_HRR$GENE %in% sub)

#####non-HRRm violin plot


all_HRD_plots_1_nonHRR<-all_HRD_plots_1[all_HRD_plots_1$Status =="Non-HRRm",]


t.test(all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="OPINION"],all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="PAOLA"])
t.test(all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="OPINION"],all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="OPINION"])
t.test(all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="OPINION"],all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="LIGHT"])
t.test(all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="OPINION"],all_HRD_plots_1_nonHRR$HRD_Score[all_HRD_plots_1_nonHRR$Study =="Study 19"])


ggplot(all_HRD_plots_1_nonHRR, aes(x = Study, y = HRD_Score),color="black")+
  geom_violin(aes(fill=Status_2),outlier.shape = NA)+
  xlab("")+
  #coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  scale_fill_manual(values=c("BRCA1m"="firebrick1","BRCA2m"="deepskyblue2",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=20, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=20),axis.title.x=element_text(size=20,face="bold"))+
  theme(axis.text.y=element_text(size=20),axis.title.y =element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))



###PAOLA
all_HRD_nonHRRm<-all_HRD[all_HRD$Status =="Non-HRRm",]
PAOLA_DTS<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/PAOLA/Copy of PAOLA_LAB_MYRIAD_20190417.csv")
PAOLA_DTS_nonHRR<-PAOLA_DTS[PAOLA_DTS$SUBJID %in% all_HRD_nonHRRm$SUBJECT[all_HRD_nonHRRm$Study =="PAOLA"],]
PAOLA_DTS_nonHRR_del<-PAOLA_DTS_nonHRR[!PAOLA_DTS_nonHRR$CLASS =="UNCERTAIN",]
PAOLA_DTS_nonHRR_del_1<-PAOLA_DTS_nonHRR_del[,c(1,14,15,7)]
colnames(PAOLA_DTS_nonHRR_del_1)<-c("SUBJECT","GENE","Mutation","HRD_Score")
PAOLA_DTS_nonHRR_del_1$Status <- "Non-HRRm"
PAOLA_DTS_nonHRR_del_1$Study <- "PAOLA"

PAOLA_DTS_nonHRR_del_1$dup<-duplicated(PAOLA_DTS_nonHRR_del_1$SUBJECT) | duplicated(PAOLA_DTS_nonHRR_del_1$SUBJECT, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique

table(PAOLA_DTS_nonHRR_del_1$dup)

PAOLA_DTS_nonHRR_del_1$dup_key<-paste(PAOLA_DTS_nonHRR_del_1$SUBJECT,PAOLA_DTS_nonHRR_del_1$GENE,sep="_")
PAOLA_DTS_nonHRR_del_1$dup_gene<-duplicated(PAOLA_DTS_nonHRR_del_1$dup_key) | duplicated(PAOLA_DTS_nonHRR_del_1$dup_key, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique
table(PAOLA_DTS_nonHRR_del_1$dup_gene)

PAOLA_DTS_nonHRR_del_1$dup_gene<-as.factor(PAOLA_DTS_nonHRR_del_1$dup_gene)
PAOLA_DTS_nonHRR_del_1$dup<-as.factor(PAOLA_DTS_nonHRR_del_1$dup)

PAOLA_DTS_nonHRR_del_1$GENE_2<-ifelse(PAOLA_DTS_nonHRR_del_1$dup =="TRUE" & PAOLA_DTS_nonHRR_del_1$dup_gene =="TRUE", as.character(PAOLA_DTS_nonHRR_del_1$GENE),
                                            ifelse(PAOLA_DTS_nonHRR_del_1$dup == "TRUE" & PAOLA_DTS_nonHRR_del_1$dup_gene =="FALSE","Co-occur with TP53",as.character(PAOLA_DTS_nonHRR_del_1$GENE)))

PAOLA_DTS_nonHRR_del_1$GENE_3<-ifelse(PAOLA_DTS_nonHRR_del_1$GENE_2 =="Co-occur with TP53","Co-occur with TP53","Unique")



###OPINION
OPINION_DTS<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/OPINION/Data_Transfer_Specs/opinion_MGSTA_randomised.csv")
OPINION_DTS_nonHRR<-OPINION_DTS[OPINION_DTS$SUBJECT %in% all_HRD_nonHRRm$SUBJECT[all_HRD_nonHRRm$Study =="OPINION"],]

del<-c(10,11)

OPINION_DTS_nonHRR_del<-OPINION_DTS_nonHRR[OPINION_DTS_nonHRR$VCLS %in% del,]
OPINION_DTS_nonHRR_del_1<-OPINION_DTS_nonHRR_del[,c(1,17,29,57)]

colnames(OPINION_DTS_nonHRR_del_1)<-c("SUBJECT","GENE","Mutation","HRD_Score")
OPINION_DTS_nonHRR_del_1$Status <- "Non-HRRm"
OPINION_DTS_nonHRR_del_1$Study <- "OPINION"


OPINION_DTS_nonHRR_del_1$dup<-duplicated(OPINION_DTS_nonHRR_del_1$SUBJECT) | duplicated(OPINION_DTS_nonHRR_del_1$SUBJECT, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique

table(OPINION_DTS_nonHRR_del_1$dup)

OPINION_DTS_nonHRR_del_1$dup_key<-paste(OPINION_DTS_nonHRR_del_1$SUBJECT,OPINION_DTS_nonHRR_del_1$GENE,sep="_")
OPINION_DTS_nonHRR_del_1$dup_gene<-duplicated(OPINION_DTS_nonHRR_del_1$dup_key) | duplicated(OPINION_DTS_nonHRR_del_1$dup_key, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique
table(OPINION_DTS_nonHRR_del_1$dup_gene)

OPINION_DTS_nonHRR_del_1$dup_gene<-as.factor(OPINION_DTS_nonHRR_del_1$dup_gene)
OPINION_DTS_nonHRR_del_1$dup<-as.factor(OPINION_DTS_nonHRR_del_1$dup)

OPINION_DTS_nonHRR_del_1$GENE_2<-ifelse(OPINION_DTS_nonHRR_del_1$dup =="TRUE" & OPINION_DTS_nonHRR_del_1$dup_gene =="TRUE", as.character(OPINION_DTS_nonHRR_del_1$GENE),
                                      ifelse(OPINION_DTS_nonHRR_del_1$dup == "TRUE" & OPINION_DTS_nonHRR_del_1$dup_gene =="FALSE","Co-occur with TP53",as.character(OPINION_DTS_nonHRR_del_1$GENE)))

OPINION_DTS_nonHRR_del_1$GENE_3<-ifelse(OPINION_DTS_nonHRR_del_1$GENE_2 =="Co-occur with TP53","Co-occur with TP53","Unique")


###LIGHT
LIGHT_DTS<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/LIGHT/light_data_consolidated_20200415.csv")
LIGHT_DTS_nonHRR<-LIGHT_DTS[LIGHT_DTS$SUBJECT %in% all_HRD_nonHRRm$SUBJECT[all_HRD_nonHRRm$Study =="LIGHT"],]
LIGHT_DTS_nonHRR_del<-LIGHT_DTS_nonHRR[LIGHT_DTS_nonHRR$VCLS %in% del,]

LIGHT_DTS_nonHRR_del_1<-LIGHT_DTS_nonHRR_del[,c(2,17,29,65)]

colnames(LIGHT_DTS_nonHRR_del_1)<-c("SUBJECT","GENE","Mutation","HRD_Score")
LIGHT_DTS_nonHRR_del_1$Status <- "Non-HRRm"
LIGHT_DTS_nonHRR_del_1$Study <- "LIGHT"

LIGHT_DTS_nonHRR_del_1$dup<-duplicated(LIGHT_DTS_nonHRR_del_1$SUBJECT) | duplicated(LIGHT_DTS_nonHRR_del_1$SUBJECT, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique

table(LIGHT_DTS_nonHRR_del_1$dup)

LIGHT_DTS_nonHRR_del_1$dup_key<-paste(LIGHT_DTS_nonHRR_del_1$SUBJECT,LIGHT_DTS_nonHRR_del_1$GENE,sep="_")
LIGHT_DTS_nonHRR_del_1$dup_gene<-duplicated(LIGHT_DTS_nonHRR_del_1$dup_key) | duplicated(LIGHT_DTS_nonHRR_del_1$dup_key, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique
table(LIGHT_DTS_nonHRR_del_1$dup_gene)

LIGHT_DTS_nonHRR_del_1$dup_gene<-as.factor(LIGHT_DTS_nonHRR_del_1$dup_gene)
LIGHT_DTS_nonHRR_del_1$dup<-as.factor(LIGHT_DTS_nonHRR_del_1$dup)

LIGHT_DTS_nonHRR_del_1$GENE_2<-ifelse(LIGHT_DTS_nonHRR_del_1$dup =="TRUE" & LIGHT_DTS_nonHRR_del_1$dup_gene =="TRUE", as.character(LIGHT_DTS_nonHRR_del_1$GENE),
                                        ifelse(LIGHT_DTS_nonHRR_del_1$dup == "TRUE" & LIGHT_DTS_nonHRR_del_1$dup_gene =="FALSE","Co-occur with TP53",as.character(LIGHT_DTS_nonHRR_del_1$GENE)))

LIGHT_DTS_nonHRR_del_1$GENE_3<-ifelse(LIGHT_DTS_nonHRR_del_1$GENE_2 =="Co-occur with TP53","Co-occur with TP53","Unique")


all_HRD_nonHRRm_groups<-rbind(PAOLA_DTS_nonHRR_del_1,OPINION_DTS_nonHRR_del_1,LIGHT_DTS_nonHRR_del_1)

write.csv(all_HRD_nonHRRm_groups,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20200825.csv",row.names=FALSE)



###study19
#Study19_DTS<-read.table("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Study19/AZO_gene_panel_DM.txt",sep="\t",header=TRUE)
Study19_nonHRR<-all_HRD_nonHRRm[all_HRD_nonHRRm$Study =="Study 19",]

Study19_nonHRR_1<-Study19_nonHRR %>% separate(SUBJECT, c("ExID1", "SUBJECT"))

study19_DTS<-read.table("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Study19/Olaparib_study19_all.data.txt",sep="\t",header=TRUE)
names(study19_DTS)[names(study19_DTS) == 'Sample.ID'] <- 'SUBJECT'
study19_DTS_nonHRR<-study19_DTS[study19_DTS$SUBJECT %in% Study19_nonHRR_1$SUBJECT,]




all_HRD_nonHRRm_19_1<-merge(Study19_nonHRR_1,study19_DTS,by = "SUBJECT",all.x=TRUE)


CLIA<-c("known","likely")
all_HRD_nonHRRm_19_1<-all_HRD_nonHRRm_19_1[all_HRD_nonHRRm_19_1$SOMATIC.STATUS.FUNCTIONAL.IMPACT %in% CLIA,]

write.csv(all_HRD_nonHRRm_19_1,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Study19/Olaparib_study19_all.data_non_HRR.csv",row.names=FALSE)

all_HRD_nonHRRm_19_2<-all_HRD_nonHRRm_19_1[,c(1,17,20,5)]

colnames(all_HRD_nonHRRm_19_2)<-c("SUBJECT","GENE","Mutation","HRD_Score")
all_HRD_nonHRRm_19_2$Status <- "Non-HRRm"
all_HRD_nonHRRm_19_2$Study <- "Study 19"

S19_DTS_nonHRR_del_1<-all_HRD_nonHRRm_19_2


S19_DTS_nonHRR_del_1$dup<-duplicated(S19_DTS_nonHRR_del_1$SUBJECT) | duplicated(S19_DTS_nonHRR_del_1$SUBJECT, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique

table(S19_DTS_nonHRR_del_1$dup)

S19_DTS_nonHRR_del_1$dup_key<-paste(S19_DTS_nonHRR_del_1$SUBJECT,S19_DTS_nonHRR_del_1$GENE,sep="_")
S19_DTS_nonHRR_del_1$dup_gene<-duplicated(S19_DTS_nonHRR_del_1$dup_key) | duplicated(S19_DTS_nonHRR_del_1$dup_key, fromLast=TRUE) ##TRUE equals co-occurring, FALSE = unique
table(S19_DTS_nonHRR_del_1$dup_gene)

S19_DTS_nonHRR_del_1$dup_gene<-as.factor(S19_DTS_nonHRR_del_1$dup_gene)
S19_DTS_nonHRR_del_1$dup<-as.factor(S19_DTS_nonHRR_del_1$dup)

S19_DTS_nonHRR_del_1$GENE_2<-ifelse(S19_DTS_nonHRR_del_1$dup =="TRUE" & S19_DTS_nonHRR_del_1$dup_gene =="TRUE", as.character(S19_DTS_nonHRR_del_1$GENE),
                                      ifelse(S19_DTS_nonHRR_del_1$dup == "TRUE" & S19_DTS_nonHRR_del_1$dup_gene =="FALSE","Co-occur with TP53",as.character(S19_DTS_nonHRR_del_1$GENE)))

S19_DTS_nonHRR_del_1$GENE_3<-ifelse(S19_DTS_nonHRR_del_1$GENE_2 =="Co-occur with TP53","Co-occur with TP53","Unique")


head(S19_DTS_nonHRR_del_1)

all_HRD_nonHRRm_groups<-rbind(PAOLA_DTS_nonHRR_del_1,OPINION_DTS_nonHRR_del_1,LIGHT_DTS_nonHRR_del_1,S19_DTS_nonHRR_del_1)

write.csv(all_HRD_nonHRRm_groups,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20200825.csv",row.names=FALSE)

gene_count<-data.frame(table(all_HRD_nonHRRm_groups$GENE))
colnames(gene_count)<-c("GENE","COUNT")
all_HRD_nonHRRm_groups<-merge(all_HRD_nonHRRm_groups,gene_count,by="GENE")

gene_count_2<-data.frame(table(all_HRD_nonHRRm_groups$GENE[all_HRD_nonHRRm_groups$HRD_Score!="FAILED"]))
colnames(gene_count_2)<-c("GENE","COUNT_withHRDscore")
all_HRD_nonHRRm_groups<-merge(all_HRD_nonHRRm_groups,gene_count_2,by="GENE")



#all_HRD_nonHRRm_groups<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20200825.csv")
all_HRD_nonHRRm_groups<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20220516.csv")

CCNE1P<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/PAOLA/PAOLA CCNE1 & ERBB2 Amplified samples.csv")

CCNE1P$nonBRCA_non_HRRm<-case_when(CCNE1P$ExId1 %in% all_HRD_nonHRRm_groups$SUBJECT~"Yes")

CCNE1P_Y<-subset(CCNE1P,CCNE1P$nonBRCA_non_HRRm =="Yes")

subset(all_HRD,all_HRD$SUBJECT %in% CCNE1P_Y$ExId1)

write.csv(CCNE1P,file= "/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/PAOLA/PAOLA CCNE1 & ERBB2 Amplified samples_nonBRCA_nonHRRm.csv",row.names=FALSE )


all_HRD_nonHRRm_groups_HRD<-all_HRD_nonHRRm_groups[!all_HRD_nonHRRm_groups$HRD_Score =="FAILED",]

all_HRD_nonHRRm_groups_HRD$Study <- factor(all_HRD_nonHRRm_groups_HRD$Study, levels = c("PAOLA","OPINION","LIGHT","Study 19"))

all_HRD_nonHRRm_groups_HRD$GENE_4 <- factor(all_HRD_nonHRRm_groups_HRD$GENE_4, levels = c("TP53","NF1","RB1","CCNE1","PIK3CA","ARID1A","KRAS","PTEN","MYH","CSMD3","MYC",
                                                                                  "MAP2K4","STAG2","ERBB2","KMT2D","TSC1","FANCM","NBN","NTHL1","PTPRD","BLM","NRAS","SOX2","MAP3K1","MCL1","Other"))
all_HRD_nonHRRm_groups_HRD$HRD_Score <-as.numeric(as.character(all_HRD_nonHRRm_groups_HRD$HRD_Score))

ggplot(all_HRD_nonHRRm_groups_HRD, aes(x = GENE_4, y = HRD_Score),color="black")+
  geom_boxplot(outlier.shape = NA)+
  xlab("GENE")+ 
  #coord_flip()+
  geom_jitter(aes(color= Study),shape=16, position=position_jitter(0.2),size=2)+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study.x,scales="free")+
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  #scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.y=element_text(size=22),axis.title.y=element_text(size=22,face="bold"))+
  theme(axis.text.x=element_text(size=12,face="bold"),axis.title.x=element_text(size=1,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')
  guides(fill=guide_legend(title=""))
  
  
  


all_HRD_nonHRRm_groups_HRD$GENE_5 <- factor(all_HRD_nonHRRm_groups_HRD$GENE_5, levels = c("TP53","Co-occur with TP53","NF1","RB1","PIK3CA","ARID1A","KRAS","PTEN","MYH","CSMD3","MYC",
                                                                                          "MAP2K4","STAG2","KMT2D","TSC1","FANCM","NBN","NTHL1","CCNE1","PTPRD","BLM","NRAS","SOX2","ERBB2","MAP3K1","MCL1","Other"))

ggplot(all_HRD_nonHRRm_groups_HRD, aes(x = GENE_5, y = HRD_Score),color="black")+
  geom_boxplot(outlier.shape = NA)+
  xlab("GENE")+ 
  #coord_flip()+
  geom_jitter(aes(color= Study),shape=16, position=position_jitter(0.2),size=2)+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Study.x,scales="free")+
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  #scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=15, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.y=element_text(size=22),axis.title.y=element_text(size=22,face="bold"))+
  theme(axis.text.x=element_text(size=12,face="bold"),axis.title.x=element_text(size=1,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')
guides(fill=guide_legend(title=""))


write.csv(all_HRD_nonHRRm_groups_HRD,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20230721.csv")


#####append all HRRm
#all_HRD<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20201020_int.csv")

#all_HRD_nonHRRm_groups<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20220516.csv")
all_HRD_nonHRRm_groups<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_nonHRRm_20230721.csv")

all_HRD_plots<-all_HRD[!all_HRD$HRD_Score =="FAILED",]
all_HRD_plots$HRD_Score <-as.numeric(as.character(all_HRD_plots$HRD_Score))
all_HRD_plots$Status <- factor(all_HRD_plots$Status, levels=c( "tBRCAm","HRRm and non-BRCAm","Non-HRRm")) #reordering zygosity variables to appear in this order in the barplot
all_HRD_plots$Status_3 <- factor(all_HRD_plots$Status_3, levels=c( "tBRCAm","Non-tBRCAm")) #reordering zygosity variables to appear in this order in the barplot

all_HRD_plots$Study <- factor(all_HRD_plots$Study, levels=c( "PAOLA","OPINION","LIGHT","Study 19", "SOLO1", "SOLO2")) #reordering zygosity variables to appear in this order in the barplot



all_HRD_plots_1_BRCAm<-all_HRD_plots[all_HRD_plots$Status =="tBRCAm",]

all_HRD_plots_1_nonHRR<-all_HRD_nonHRRm_groups


all_HRD_plots_1_nonHRR<-all_HRD_plots[all_HRD_plots$Status =="Non-HRRm",]

all_HRD_plots_1_nonHRR$key<- paste(all_HRD_plots_1_nonHRR$SUBJECT,all_HRD_plots_1_nonHRR$Study,sep="")
all_HRD_plots_1_BRCAm$key<- paste(all_HRD_plots_1_BRCAm$SUBJECT,all_HRD_plots_1_BRCAm$Study,sep="")

all_HRD_plots_1_nonHRR<-all_HRD_plots_1_nonHRR[!duplicated(all_HRD_plots_1_nonHRR$key),]
all_HRD_plots_1_BRCAm<-all_HRD_plots_1_BRCAm[!duplicated(all_HRD_plots_1_BRCAm$key),]

all_HRD_plots_1_nonHRR$GENE_4<-case_when(is.na(all_HRD_plots_1_nonHRR$GENE) ~"non-HRRm")
all_HRD_plots_1_BRCAm$GENE_4<-case_when(all_HRD_plots_1_BRCAm$Status =="tBRCAm" ~"tBRCAm")


shared_cols<-intersect(colnames(all_HRD_plots_1_nonHRR),colnames(all_HRD_nonHRRm_groups_HRD))

all_HRD_plots_1_nonHRR_1<-all_HRD_plots_1_nonHRR[ , which(colnames(all_HRD_plots_1_nonHRR) %in% shared_cols)]
all_HRD_plots_1_BRCAm_1<-all_HRD_plots_1_BRCAm[ , which(colnames(all_HRD_plots_1_BRCAm) %in% shared_cols)]

all_HRD_nonHRRm_groups_HRD_1<-all_HRD_nonHRRm_groups_HRD[,which(colnames(all_HRD_nonHRRm_groups_HRD) %in% shared_cols)]
all_HRD_plots_1_BRCAm_1$Status_2<-"tBRCAm"
all_HRD_plots_1_nonHRR_1$Status_2<-"non-HRRm"
all_HRD_nonHRRm_groups_HRD_1$Status_2 <- "non-HRRm (Gene by Gene)"

all_HRD_plots_1_BRCAm_1<-all_HRD_plots_1_BRCAm_1[!is.na(all_HRD_plots_1_BRCAm_1$SUBJECT),]

B_1<-rbind(all_HRD_plots_1_BRCAm_1,all_HRD_plots_1_nonHRR_1,all_HRD_nonHRRm_groups_HRD_1)

B_1$Study <- factor(B_1$Study, levels = c("PAOLA","OPINION","LIGHT","Study 19"))
B_1$Status_2 <- gsub("non-HRRm", "Non-HRRm",B_1$Status_2)
B_1$Status_2 <- gsub("non-HRRm (Gene by Gene)", "Non-HRRm (Gene by Gene)",B_1$Status_2)

B_1$Status_2 <- factor(B_1$Status_2, levels = c("tBRCAm","Non-HRRm","Non-HRRm (Gene by Gene)"))

#B_1$GENE_4 <- factor(B_1$GENE_4, levels = c("tBRCAm","non-HRRm","TP53","NF1","RB1","CCNE1","PIK3CA","ARID1A","KRAS","PTEN","MYH","CSMD3","MYC",
                                                                                          "MAP2K4","STAG2","ERBB2","KMT2D","TSC1","FANCM","NBN","NTHL1","PTPRD","BLM","NRAS","SOX2","MAP3K1","MCL1","Other"))

B_1$GENE_4 <- gsub("non-HRRm", "Non-HRRm",B_1$GENE_4)

B_1$GENE_4 <- factor(B_1$GENE_4, levels = c("tBRCAm","Non-HRRm","MAP3K1","MYC","RB1","KMT2D","CSMD3","SOX2","NF1","MAP2K4","NTHL1","PTEN","STAG2","PTPRD","TP53","FANCM","ERBB2","TSC1","MCL1","BLM","PIK3CA","MYH","NBN","CCNE1","ARID1A","NRAS","KRAS","Other"))


B_1$HRD_Score <-as.numeric(as.character(B_1$HRD_Score))

B_1<-B_1[!is.na(B_1$Study),]


  
ggplot(B_1, aes(x = GENE_4, y = HRD_Score),color="black")+
#ggplot(B_1, aes(x = fct_reorder(GENE_4, HRD_Score, .desc =TRUE), y= HRD_Score))+
  geom_boxplot(outlier.shape = NA)+
  xlab("GENE")+ 
  #coord_flip()+
  geom_jitter(aes(color= Study),shape=16, position=position_jitter(0.2),size=2)+
  geom_hline(yintercept=42,linetype="dashed")+
  facet_grid(.~ Status_2,scales="free",space="free")+
  theme(panel.spacing = unit(4, "lines"))+
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  #scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=13, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.y=element_text(size=22),axis.title.y=element_text(size=22,face="bold"))+
  theme(axis.text.x=element_text(size=12,face="bold"),axis.title.x=element_text(size=1,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')
  guides(fill=guide_legend(title=""))

ggplot(B_1, aes(x = GENE_4, y = HRD_Score),color="black")+
  #ggplot(B_1, aes(x = fct_reorder(GENE_4, HRD_Score, .desc =TRUE), y= HRD_Score))+
  geom_boxplot(outlier.shape = NA,aes(fill=GENE_4))+
  scale_fill_manual(values=as.vector(polychrome(28)))+
    xlab("GENE")+ 
  #coord_flip()+
  geom_jitter(color= "black",shape=16, position=position_jitter(0.2),size=1)+
  geom_hline(yintercept=42,linetype="dashed")+
  facet_grid(.~ Status_2,scales="free",space="free")+
  
  scale_y_continuous(limits = c(-1, 102), breaks = seq(0, 103, by = 25))+
  #scale_color_manual(values=c("Biallelic Inactivation"="blue3", "Heterozygous"="grey","Unknown"="grey5"))+
  theme_classic()+
  theme(panel.spacing = unit(1, "lines"))+
  
  theme(strip.text.x = element_text(size=12, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.y=element_text(size=22),axis.title.y=element_text(size=22,face="bold"))+
  theme(axis.text.x=element_text(size=20,face="bold",angle =45,,hjust=1),axis.title.x=element_blank())+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("GIS")+
  theme(legend.position = 'none')
guides(fill=guide_legend(title=""))





####pie chart


test2<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/test_PAOLA_OPINION_LIGHT_S19.csv")


percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC((x/1788)*100, format = format, digits = digits, ...), "%")
}


test2$fill2 <- factor(test2$fill2, levels = c(
  "HRD negative",
  "HRD negative (HRRm)",
  "HRD negative (non-HRRm)",
  "HRD unknown",
  "HRD unknown (HRRm)",
  "HRD unknown (non-HRRm)",
  "HRD positive",
  "tBRCAm",
  "HRD positive (HRRm)",
  "HRD positive (non-HRRm)"
))

tiff("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study//pie_HRD-HRR_POLS19_blue.tiff", units="in", width=50, height=50, res=100)
ggplot(test2, aes(x = level, y = Value, fill = fill2))+
  geom_col(width = 1, color = "white", size = 0.25, position = position_stack())+
  coord_polar(theta = "y")+
  #geom_text(aes(label = Value), size = 6,color="white", position = position_stack(vjust = 0.5))+
  #geom_text(aes(label = paste0(Value,"\n", percent(Value))), size = 5,color="white", position = position_stack(vjust = 0.45))+
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL)+
  scale_fill_manual(values=c("HRD negative"="snow3",
                             "HRD negative (HRRm)" = "orange",
                             "HRD negative (non-HRRm)" ="snow3",
                             "HRD positive" ="royalblue3",
                             "tBRCAm"="royalblue4",
                             "HRD positive (HRRm)"="orange",
                             "HRD positive (non-HRRm)"="royalblue1",
                             "HRD unknown"="snow4",
                             "HRD unknown (HRRm)" = "orange",
                             "HRD unknown (non-HRRm)" = "snow4"), na.translate = F)+
  #scale_alpha_manual(values = c( "1" = 1, "2" = 0.55), guide = F)
  labs(x = NULL, y = NULL)+ 
  theme_void()+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title="HRD status"))+
  theme(legend.title=element_text(size=12,face="bold"))+
  theme(legend.text=element_text(size=12))
dev.off()




##pie chart HRR
prev4<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/HRRm_non_BRCAm_transposed.csv")

#prev4<-prev4[!prev4$Test =="PROpel_Tissue",]
mycols2<-c("navy","forestgreen","purple","lightcoral","deeppink","dodgerblue2","brown2","cyan","darkorange1","blue","darkkhaki","yellow","darkred","darkgoldenrod1","turquoise4","maroon4","orange")
prev4$Study <- factor(prev4$Study, levels = c("PAOLA-1 (n=54)","OPINION (n=34)","LIGHT (n=21)","ORZORA (n=33)","Study 19 (n=21)"))

prev4$GENE <- factor(prev4$GENE, levels = c("Co-occurring genes",
                                            "FANCA only",
                                            "FANC2D only",
                                            "RAD52 only",
                                            "XRCC3 only",
                                            "FANCL only",
                                            "RAD51D only",
                                            "BARD1 only",
                                            "RAD51B only",
                                            "RAD51C only",
                                            "RAD54L only",
                                            "PALB2 only",
                                            "BRIP1 only",
                                            "PPP2R2A only",
                                            "CHEK2 only",
                                            "CDK12 only",
                                            "ATM only"))


tiff("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Figures/pie_HRRm_non_BRCA_cross_study.tiff", units="in", width=50, height=50, res=100)
ggplot(prev4, aes(x = factor(1), y = Percent, fill = GENE)) +
  geom_bar(width = 1, size=1,stat = "identity", color = "black",alpha=0.8)+
  facet_grid(.~ Study)+
  coord_polar(theta="y")+
  theme(axis.text.x=element_text(color='black'))+
  #geom_text(aes(label = paste0(Count,"\n (",  Frequency....,"%",")")), size = 4,color="white", position = position_stack(vjust = 0.45))
  scale_fill_manual(values = mycols2,breaks=c("ATM only","CDK12 only","CHEK2 only","BRIP1 only","PALB2 only","RAD54L only","RAD51B only","RAD51C only","BARD1 only","RAD51D only","FANCL only","FANCA only",
                                              "FANC2D only",
                                              "RAD52 only",
                                              "XRCC3 only", "PPP2R2A only","Co-occurring genes"))+
  theme_void()+
  theme(strip.text.x = element_text(size=22, color="Black",
                                    face="bold.italic"))+
  theme(panel.spacing = unit(2, "lines"))+
  guides(fill=guide_legend(title="HRR genes"))+
  #theme(strip.background = element_rect(colour="black", fill="white", size=2, linetype="solid"))+
  #theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+
  theme(legend.title=element_text(size=20,face="bold"))+
  theme(legend.text=element_text(size=20))+
  #scale_y_continuous(limits = c(0,101))+
  theme(legend.position = 'right')
dev.off()


