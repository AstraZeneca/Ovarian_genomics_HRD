library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggpubr")
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
#install.packages("wesanderson")
library(wesanderson)
library(viridis)
library(pals)

all_HRD<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_20201020_int_base.csv")



###baseline_char
b_s19<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_S19.csv")
s19_int<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Study19/Study19_new.csv")

b_PAOLA<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_PAOLA1.csv")
b_S1<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_SOLO1.csv")
b_S2<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_SOLO2.csv")
h_S1<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_rsfa_SOLO1.csv")
h_S2<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_rsfa_SOLO2.csv")
h_L2<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_admh_LIGHT.csv")
b_OP<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_OPINION.csv")
b_L<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/Baseline_Characteristics/deid_adsl_LIGHT.csv")




#RACE

####study 19
b_s19_sub<-b_s19[!b_s19$ARM =="Screen Failure",]
dim(b_s19_sub)

all_HRD_s19<-subset(all_HRD,all_HRD$Study =="Study 19")
all_HRD_s191<-all_HRD_s19 %>% separate(SUBJECT, c("Specimen_ID", "ECODE"))

int_deid_HRD_s19<-intersect(all_HRD_s191$ECODE,b_s19_sub$SUBJID)

test<-setdiff(all_HRD_s191$ECODE,b_s19_sub$SUBJID) 

s19_int_excl_deid<-subset(s19_int,s19_int$Sample %in% test) 

table(s19_int_excl_deid$COUNTRY)




####

data.frame(table(b_s19_sub$RACE))

data.frame(table(b_s19_sub$RACE[b_s19_sub$SUBJID %in% int_deid_HRD_s19]))


b_s19_sub$RACE_AB<-ifelse(b_s19_sub$RACE =="WHITE", "White",
                         ifelse(b_s19_sub$RACE =="ASIAN", "Asian",
                                ifelse(b_s19_sub$RACE =="BLACK OR AFRICAN AMERICAN",  "Black or African American","Other/Unknown")))

data.frame(table(b_s19_sub$RACE_AB[b_s19_sub$SUBJID %in% int_deid_HRD_s19]))


####HISTOLOGY

b_s19_sub$HIST_AB<-"Serous"



data.frame(table(b_s19_sub$LOCATION[b_s19_sub$SUBJID %in% int_deid_HRD_s19])) 

b_s19_sub$LOCATION_AB<-ifelse(b_s19_sub$LOCATION =="OVARY", "Ovary",
                          ifelse(b_s19_sub$LOCATION =="FALLOPIAN TUBE", "Fallopian Tube or Primary Peritoneal",
                                 ifelse(b_s19_sub$LOCATION =="PRIMARY PERITONEAL", "Fallopian Tube or Primary Peritoneal","Other/Unknown")))

b_s19_sub$LOCATION_AB<-ifelse(b_s19_sub$LOCATION =="OVARY", "Ovary",
                              ifelse(b_s19_sub$LOCATION =="FALLOPIAN TUBE", "Fallopian Tube",
                                     ifelse(b_s19_sub$LOCATION =="PRIMARY PERITONEAL", "Peritoneal","Other/Unknown")))


data.frame(table(b_s19_sub$LOCATION_AB[b_s19_sub$SUBJID %in% int_deid_HRD_s19]))



data.frame(table(b_s19_sub$FIGOSTG[b_s19_sub$SUBJID %in% int_deid_HRD_s19])) ###FIGO

b_s19_sub$FIGO_AB<-ifelse(b_s19_sub$FIGOSTG =="III", "III",
                          ifelse(b_s19_sub$FIGOSTG =="II", "II",
                                 ifelse(b_s19_sub$FIGOSTG =="IIB", "II",
                                        ifelse(b_s19_sub$FIGOSTG =="IIC", "II",
                                               ifelse(b_s19_sub$FIGOSTG =="IIIA", "III",
                                                      ifelse(b_s19_sub$FIGOSTG =="IIIB", "III",
                                                             ifelse(b_s19_sub$FIGOSTG =="IIIC", "III",
                                                                    ifelse(b_s19_sub$FIGOSTG =="IV", "IV", "Other/Unknown"))))))))


data.frame(table(b_s19_sub$FIGO_AB[b_s19_sub$SUBJID %in% int_deid_HRD_s19])) ###FIGO

b_s19_sub_sub<-b_s19_sub[,c(4,83:86)]
b_s19_sub_sub$Study<-"Study 19"

####PAOLA-1 
####RACE cannot be done as this information was not collected
b_PAOLA_sub<-b_PAOLA[!b_PAOLA$ARM =="SCREEN FAILURE",]
b_PAOLA_sub$RACE_AB<-"NA"

####HISTOLOGY

dim(b_PAOLA_sub)


all_HRD_P<-subset(all_HRD,all_HRD$Study =="PAOLA")

length(intersect(all_HRD_P$SUBJECT, b_PAOLA_sub$SUBJID))


data.frame(table(b_PAOLA_sub$HSTLTYP)) ###Histology



b_PAOLA_sub$HIST_AB<-ifelse(b_PAOLA_sub$HSTLTYP =="SEROUS", "Serous",
                              ifelse(b_PAOLA_sub$HSTLTYP =="ENDOMETRIOID", "Endometrioid","Other/Unknown"))

data.frame(table(b_PAOLA_sub$HIST_AB)) ###Histology


data.frame(table(b_PAOLA_sub$PTUMLOC)) ###Location



#b_PAOLA_sub$LOCATION_AB<-ifelse(b_PAOLA_sub$PTUMLOC =="Ovary", "Ovary",
 #                             ifelse(b_PAOLA_sub$PTUMLOC =="Fallopian", "Fallopian Tube or Primary Peritoneal",
  #                                   ifelse(b_PAOLA_sub$PTUMLOC =="Peritoneal", "Fallopian Tube or Primary Peritoneal","Other/Unknown")))

b_PAOLA_sub$LOCATION_AB<-ifelse(b_PAOLA_sub$PTUMLOC =="Ovary", "Ovary",
                                ifelse(b_PAOLA_sub$PTUMLOC =="Fallopian", "Fallopian Tube",
                                       ifelse(b_PAOLA_sub$PTUMLOC =="Peritoneal", "Peritoneal","Other/Unknown")))

data.frame(table(b_PAOLA_sub$LOCATION_AB)) ###Location


data.frame(table(b_PAOLA_sub$FIGOBL)) ###FIGO




b_PAOLA_sub$FIGO_AB<-ifelse(b_PAOLA_sub$FIGOBL =="III B", "III",
                                 ifelse(b_PAOLA_sub$FIGOBL =="III C", "III",
                                        ifelse(b_PAOLA_sub$FIGOBL =="IV", "IV", "Other/Unknown")))


data.frame(table(b_PAOLA_sub$FIGO_AB))

b_PAOLA_sub_sub<-b_PAOLA_sub[,c(19,112:115)]
b_PAOLA_sub_sub$Study<-"PAOLA"



####OPINION

##RACE
b_OP_sub<-b_OP[!b_OP$ARM =="Screen Failure",]
dim(b_OP_sub)

all_HRD_O<-subset(all_HRD,all_HRD$Study =="OPINION")

length(intersect(all_HRD_O$SUBJECT,b_OP_sub$SUBJID))

data.frame(table(b_OP_sub$RACE))


data.frame(table(b_OP_sub$RACE[b_OP_sub$SU %in% all_HRD_O$SUBJECT]))


b_OP_sub$RACE_AB<-ifelse(b_OP_sub$RACE =="WHITE", "White",
                         ifelse(b_OP_sub$RACE =="ASIAN", "Asian",
                                ifelse(b_OP_sub$RACE =="BLACK OR AFRICAN AMERICAN",  "Black or African American","Other/Unknown")))

data.frame(table(b_OP_sub$RACE_AB[b_OP_sub$SUBJID %in% all_HRD_O$SUBJECT]))

###HISTOLOGY
data.frame(table(b_OP_sub$HSTYPBL)) ###Histology



b_OP_sub$HIST_AB<-ifelse(b_OP_sub$HSTYPBL =="SEROUS", "Serous",
                            ifelse(b_OP_sub$HSTYPBL =="ENDOMETRIOID", "Endometrioid","Other/Unknown"))

data.frame(table(b_OP_sub$HIST_AB)) ###Histology





data.frame(table(b_OP_sub$PTULOCBL)) ####Location



b_OP_sub$LOCATION_AB<-ifelse(b_OP_sub$PTULOCBL =="OVARY", "Ovary",
                             ifelse(b_OP_sub$PTULOCBL =="PERITONEUM", "Peritoneal",
                                    ifelse(b_OP_sub$PTULOCBL =="FALLOPIAN TUBE", "Fallopian Tube","Other/Unknown")))

data.frame(table(b_OP_sub$LOCATION_AB)) ###Location



data.frame(table(b_OP_sub$FIGOSTBL)) ####FIGO stage


b_OP_sub$FIGO_AB<-ifelse(b_OP_sub$FIGOSTBL =="STAGE III", "III",
                          ifelse(b_OP_sub$FIGOSTBL =="STAGE II", "II",
                                 ifelse(b_OP_sub$FIGOSTBL =="STAGE IIA", "II",
                                 ifelse(b_OP_sub$FIGOSTBL =="STAGE IIB", "II",
                                        ifelse(b_OP_sub$FIGOSTBL =="STAGE IIC", "II",
                                               ifelse(b_OP_sub$FIGOSTBL =="STAGE IIIA", "III",
                                                      ifelse(b_OP_sub$FIGOSTBL =="STAGE IIIB", "III",
                                                             ifelse(b_OP_sub$FIGOSTBL =="STAGE IIIC", "III",
                                                                    ifelse(b_OP_sub$FIGOSTBL =="STAGE IV", "IV",
                                                                           ifelse(b_OP_sub$FIGOSTBL =="STAGE IVA", "IV", "Other/Unknown"))))))))))


data.frame(table(b_OP_sub$FIGO_AB)) ###FIGO

b_OP_sub_sub<-b_OP_sub[,c(3,94:97)]
b_OP_sub_sub$Study<-"OPINION"


####SOLO1
##RACE
b_S1_sub<-b_S1[!b_S1$ARM =="Screen Failure",]
b_S1_sub<-b_S1_sub[!b_S1_sub$ARM =="Not Assigned",]

dim(b_S1_sub)
table(b_S1_sub$COUNTRY)
dim(b_S1_sub)
#[1] 301  62

all_HRD_S1<-subset(all_HRD,all_HRD$Study =="SOLO1")
dim(all_HRD_S1)
#[1] 284  17
length(intersect(all_HRD_S1$SUBJECT,b_S1_sub$SUBJID))

deid_bHRD_SOLO1_diff<-setdiff(all_HRD_S1$SUBJECT,b_S1_sub$SUBJID)

subset(all_HRD_S1,all_HRD_S1$SUBJECT %in% deid_bHRD_SOLO1_diff)

data.frame(table(b_S1_sub$RACE[b_S1_sub$SUBJID %in% all_HRD_S1$SUBJECT]))



b_S1_sub$RACE_AB<-ifelse(b_S1_sub$RACE =="WHITE", "White",
                          ifelse(b_S1_sub$RACE =="ASIAN", "Asian",
                                 ifelse(b_S1_sub$RACE =="BLACK OR AFRICAN AMERICAN",  "Black or African American","Other/Unknown")))

data.frame(table(b_S1_sub$RACE_AB[b_S1_sub$SUBJID %in% all_HRD_S1$SUBJECT]))

b_S1_sub_sub<-b_S1_sub[,c(4,63)]


###HISTOLOGY
####Can't find tumour histology in ADSL file but FIGO stage is present. Information is in primary manuscript though so information must have been collected
###histology is rsfa file

h_S1_sub<-h_S1[!h_S1$ARM =="Screen Failure",]
length(intersect(h_S1_sub$SUBJID,all_HRD_S1$SUBJECT))

h_S1_sub_hist<-subset(h_S1_sub,h_S1_sub$FATEST == "Histology Type")
dim(h_S1_sub_hist)

data.frame(table(h_S1_sub_hist$AVALC[h_S1_sub_hist$SUBJID %in% all_HRD_S1$SUBJECT]))

h_S1_sub_hist$HIST_AB<-ifelse(h_S1_sub_hist$AVALC =="SEROUS", "Serous",
                         ifelse(h_S1_sub_hist$AVALC =="ENDOMETROID", "Endometrioid","Other/Unknown"))

data.frame(table(h_S1_sub_hist$HIST_AB[h_S1_sub_hist$SUBJID %in% all_HRD_S1$SUBJECT])) ###Histology

h_S1_sub_hist_sub<-h_S1_sub_hist[,c(4,60)]
dim(h_S1_sub_hist_sub)

####location

h_S1_sub_loc<-subset(h_S1_sub,h_S1_sub$FATEST == "Primary Tumor Location")
dim(h_S1_sub_loc)

data.frame(table(h_S1_sub_loc$AVALC[h_S1_sub_loc$SUBJID %in% all_HRD_S1$SUBJECT]))

#h_S1_sub_loc$LOCATION_AB<-ifelse(h_S1_sub_loc$AVALC =="OVARY", "Ovary",
 #                            ifelse(h_S1_sub_loc$AVALC =="PRIMARY PERITONEAL", "Fallopian Tube or Primary Peritoneal",
  #                                  ifelse(h_S1_sub_loc$AVALC =="FALLOPIAN TUBES", "Fallopian Tube or Primary Peritoneal","Other/Unknown")))


h_S1_sub_loc$LOCATION_AB<-ifelse(h_S1_sub_loc$AVALC =="OVARY", "Ovary",
                                 ifelse(h_S1_sub_loc$AVALC =="PRIMARY PERITONEAL", "Peritoneal",
                                        ifelse(h_S1_sub_loc$AVALC =="FALLOPIAN TUBES", "Fallopian Tube","Other/Unknown")))


data.frame(table(h_S1_sub_loc$LOCATION_AB[h_S1_sub_loc$SUBJID %in% all_HRD_S1$SUBJECT])) ###Histology

h_S1_sub_loc_sub<-h_S1_sub_loc[,c(4,60)]

dim(h_S1_sub_loc_sub)

###FIGO
h_S1_sub_fgo<-subset(h_S1_sub,h_S1_sub$FATEST == "FIGO Stage")
dim(h_S1_sub_fgo)

data.frame(table(h_S1_sub_fgo$AVALC[h_S1_sub_fgo$SUBJID %in% all_HRD_S1$SUBJECT]))


h_S1_sub_fgo$FIGO_AB<-ifelse(h_S1_sub_fgo$AVALC =="Stage II", "II",
                             ifelse(h_S1_sub_fgo$AVALC =="Stage IIA", "II",
                                    ifelse(h_S1_sub_fgo$AVALC =="Stage IIB", "II",
                                           ifelse(h_S1_sub_fgo$AVALC =="Stage IIC", "II",
                                                  ifelse(h_S1_sub_fgo$AVALC =="Stage III", "III",
                                                         ifelse(h_S1_sub_fgo$AVALC =="Stage IIIA", "III",
                                                                ifelse(h_S1_sub_fgo$AVALC =="Stage IIIB", "III",
                                                                       ifelse(h_S1_sub_fgo$AVALC =="Stage IIIC", "III",
                                                                              ifelse(h_S1_sub_fgo$AVALC =="Stage IV", "IV", "Other/Unknown")))))))))


data.frame(table(h_S1_sub_fgo$FIGO_AB[h_S1_sub_fgo$SUBJID %in% all_HRD_S1$SUBJECT])) ###Histology

h_S1_sub_fgo_sub<-h_S1_sub_fgo[,c(4,60)]

dim(h_S1_sub_fgo_sub)



###need to merge all the files together (Race, Histology, Location,FIGO, Study)

b_S1_sub_sub_merge1<-merge(b_S1_sub_sub,h_S1_sub_hist_sub,by="SUBJID",all.x=TRUE)
b_S1_sub_sub_merge2<-merge(b_S1_sub_sub_merge1,h_S1_sub_loc_sub,by="SUBJID",all.x=TRUE)
b_S1_sub_sub_merge3<-merge(b_S1_sub_sub_merge2,h_S1_sub_fgo_sub,by="SUBJID",all.x=TRUE)


data.frame(table(b_S1_sub_sub_merge3$RACE_AB[b_S1_sub_sub_merge3$SUBJID %in% all_HRD_S1$SUBJECT]))
data.frame(table(b_S1_sub_sub_merge3$HIST_AB[b_S1_sub_sub_merge3$SUBJID %in% all_HRD_S1$SUBJECT]))
data.frame(table(b_S1_sub_sub_merge3$LOCATION_AB[b_S1_sub_sub_merge3$SUBJID %in% all_HRD_S1$SUBJECT]))
data.frame(table(b_S1_sub_sub_merge3$FIGO_AB[b_S1_sub_sub_merge3$SUBJID %in% all_HRD_S1$SUBJECT]))


b_S1_sub_sub_merge3$Study<-"SOLO1"


######SOLO2

###RACE
b_S2_sub<-b_S2[!b_S2$ARM =="Screen Failure",]
b_S2_sub<-b_S2_sub[!b_S2_sub$ARM =="Not Assigned",]

dim(b_S2_sub)
#[1] 228  69
table(b_S2_sub$COUNTRY)


all_HRD_S2<-subset(all_HRD,all_HRD$Study =="SOLO2")

length(intersect(all_HRD_S2$SUBJECT,b_S2_sub$SUBJID))

data.frame(table(b_S2_sub$RACE))

data.frame(table(b_S2_sub$RACE[b_S2_sub$SUBJID %in% all_HRD_S2$SUBJECT]))



b_S2_sub$RACE_AB<-ifelse(b_S2_sub$RACE =="WHITE", "White",
                         ifelse(b_S2_sub$RACE =="ASIAN", "Asian",
                                ifelse(b_S2_sub$RACE =="BLACK OR AFRICAN AMERICAN",  "Black or African American","Other/Unknown")))

data.frame(table(b_S2_sub$RACE_AB[b_S2_sub$SUBJID %in% all_HRD_S2$SUBJECT]))
b_S2_sub_sub<-b_S2_sub[,c(4,70)]

dim(b_S2_sub_sub)
#[1] 228   2


###HISTOLOGY
###histology is rsfa file

h_S2_sub<-h_S2[!h_S2$ARM =="Screen Failure",]
dim(h_S2_sub)

length(intersect(h_S2_sub$SUBJID,all_HRD_S2$SUBJECT))

h_S2_sub_hist<-subset(h_S2_sub,h_S2_sub$FATEST == "Histology Type")
dim(h_S2_sub_hist)


data.frame(table(h_S2_sub_hist$AVALC[h_S2_sub_hist$SUBJID %in% all_HRD_S2$SUBJECT]))

h_S2_sub_hist$HIST_AB<-ifelse(h_S2_sub_hist$AVALC =="SEROUS", "Serous",
                              ifelse(h_S2_sub_hist$AVALC =="ENDOMETROID", "Endometrioid","Other/Unknown"))

data.frame(table(h_S2_sub_hist$HIST_AB[h_S2_sub_hist$SUBJID %in% all_HRD_S2$SUBJECT])) ###Histology

h_S2_sub_hist_sub<-h_S2_sub_hist[,c(4,62)]
dim(h_S2_sub_hist_sub)


####location

h_S2_sub_loc<-subset(h_S2_sub,h_S2_sub$FATEST == "Primary Tumor Location")
dim(h_S2_sub_loc)

data.frame(table(h_S2_sub_loc$AVALC[h_S2_sub_loc$SUBJID %in% all_HRD_S2$SUBJECT]))

h_S2_sub_loc$LOCATION_AB<-ifelse(h_S2_sub_loc$AVALC =="OVARY", "Ovary",
                                 ifelse(h_S2_sub_loc$AVALC =="PRIMARY PERITONEAL", "Peritoneal",
                                        ifelse(h_S2_sub_loc$AVALC =="FALLOPIAN TUBES", "Fallopian Tube","Other/Unknown")))

data.frame(table(h_S2_sub_loc$LOCATION_AB[h_S2_sub_loc$SUBJID %in% all_HRD_S2$SUBJECT])) ###Histology

h_S2_sub_loc_sub<-h_S2_sub_loc[,c(4,62)]
dim(h_S2_sub_loc_sub)


###FIGO

h_S2_sub_fgo<-subset(h_S2_sub,h_S2_sub$FATEST == "FIGO Stage")
dim(h_S2_sub_fgo)

data.frame(table(h_S2_sub_fgo$AVALC[h_S2_sub_fgo$SUBJID %in% all_HRD_S2$SUBJECT]))


h_S2_sub_fgo$FIGO_AB<-ifelse(h_S2_sub_fgo$AVALC =="Stage II", "II",
                             ifelse(h_S2_sub_fgo$AVALC =="Stage IIA", "II",
                                    ifelse(h_S2_sub_fgo$AVALC =="Stage IIB", "II",
                                           ifelse(h_S2_sub_fgo$AVALC =="Stage IIC", "II",
                                                  ifelse(h_S2_sub_fgo$AVALC =="Stage III", "III",
                                                         ifelse(h_S2_sub_fgo$AVALC =="Stage IIIA", "III",
                                                                ifelse(h_S2_sub_fgo$AVALC =="Stage IIIB", "III",
                                                                       ifelse(h_S2_sub_fgo$AVALC =="Stage IIIC", "III",
                                                                              ifelse(h_S2_sub_fgo$AVALC =="Stage IV", "IV", "Other/Unknown")))))))))


data.frame(table(h_S2_sub_fgo$FIGO_AB[h_S2_sub_fgo$SUBJID %in% all_HRD_S2$SUBJECT])) ###Histology

h_S2_sub_fgo_sub<-h_S2_sub_fgo[,c(4,62)]

dim(h_S2_sub_fgo_sub)


###need to merge all the files together (Race, Histology, Location,FIGO, Study)

b_S2_sub_sub_merge1<-merge(b_S2_sub_sub,h_S2_sub_hist_sub,by="SUBJID",all.x=TRUE)
b_S2_sub_sub_merge2<-merge(b_S2_sub_sub_merge1,h_S2_sub_loc_sub,by="SUBJID",all.x=TRUE)
b_S2_sub_sub_merge3<-merge(b_S2_sub_sub_merge2,h_S2_sub_fgo_sub,by="SUBJID",all.x=TRUE)


data.frame(table(b_S2_sub_sub_merge3$RACE_AB[b_S2_sub_sub_merge3$SUBJID %in% all_HRD_S2$SUBJECT]))
data.frame(table(b_S2_sub_sub_merge3$HIST_AB[b_S2_sub_sub_merge3$SUBJID %in% all_HRD_S2$SUBJECT]))
data.frame(table(b_S2_sub_sub_merge3$LOCATION_AB[b_S2_sub_sub_merge3$SUBJID %in% all_HRD_S2$SUBJECT]))
data.frame(table(b_S2_sub_sub_merge3$FIGO_AB[b_S2_sub_sub_merge3$SUBJID %in% all_HRD_S2$SUBJECT]))


b_S2_sub_sub_merge3$Study<-"SOLO2"


b_all_sub_sub<-rbind(b_s19_sub_sub,b_PAOLA_sub_sub,b_OP_sub_sub,b_S1_sub_sub_merge3,b_S2_sub_sub_merge3)


######LIGHT
b_L_sub<-b_L[!b_L$ARM =="Screen Failure",]
dim(b_L_sub)
table(b_L_sub$COUNTRY)


all_HRD_L<-subset(all_HRD,all_HRD$Study =="LIGHT")

length(intersect(all_HRD_L$SUBJECT,b_L_sub$SUBJID))
length(intersect(all_HRD_L$SUBJECT,b_L$SUBJID))



data.frame(table(b_L_sub$RACE))

data.frame(table(b_L_sub$RACE[b_L_sub$SUBJID %in% all_HRD_L$SUBJECT]))

data.frame(table(b_L$RACE[b_L$SUBJID %in% all_HRD_L$SUBJECT]))


b_L_sub$RACE_AB<-ifelse(b_L_sub$RACE =="WHITE", "White",
                         ifelse(b_L_sub$RACE =="ASIAN", "Asian",
                                ifelse(b_L_sub$RACE =="BLACK OR AFRICAN AMERICAN",  "Black or African American","Other/Unknown")))

data.frame(table(b_L_sub$RACE_AB))
b_L_sub_sub<-b_L_sub[,c(4,80)]


####Histology
h_L2_sub<-h_L2[!h_L2$ARM =="Screen Failure",]
dim(h_L2_sub)

length(intersect(h_L2_sub$SUBJID,all_HRD_L$SUBJECT))

h_L2_sub_hist<-subset(h_L2_sub,h_L2_sub$MHTERM == "OVARIAN/FALLOPIAN/PERITONEAL CANCER")
dim(h_L2_sub_hist)
#[1] 263  47 )


data.frame(table(h_L2_sub_hist$HISTTY[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT]))

h_L2_sub_hist$HIST_AB<-ifelse(h_L2_sub_hist$HISTTYP =="SEROUS", "Serous",
                              ifelse(h_L2_sub_hist$HISTTYP =="ENDOMETRIOID", "Endometrioid","Other/Unknown"))

data.frame(table(h_L2_sub_hist$HIST_AB[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT])) ###Histology




####location


data.frame(table(h_L2_sub_hist$TUMLOC[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT]))

h_L2_sub_hist$LOCATION_AB<-ifelse(h_L2_sub_hist$TUMLOC =="OVARY", "Ovary",
                                 ifelse(h_L2_sub_hist$TUMLOC =="PRIMARY PERITONEAL", "Peritoneal",
                                        ifelse(h_L2_sub_hist$TUMLOC =="FALLOPIAN TUBE", "Fallopian Tube","Other/Unknown")))

data.frame(table(h_L2_sub_hist$LOCATION_AB[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT])) ###Histology




###FIGO


data.frame(table(h_L2_sub_hist$FIGOSTAG[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT]))


h_L2_sub_hist$FIGO_AB<-ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE II", "II",
                             ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIA", "II",
                                    ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIB", "II",
                                           ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIC", "II",
                                                  ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE III", "III",
                                                         ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIIA", "III",
                                                                ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIIB", "III",
                                                                       ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IIIC", "III",
                                                                              ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IV", "IV",
                                                                                     ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IVA", "IV",
                                                                                            ifelse(h_L2_sub_hist$FIGOSTAG =="STAGE IVB","IV","Other/Unknown")))))))))))



data.frame(table(h_L2_sub_hist$FIGO_AB[h_L2_sub_hist$SUBJID %in% all_HRD_L$SUBJECT])) ###Histology


h_L_sub_sub<- h_L2_sub_hist[,c(4,48:50)]
b_L_sub_sub1<-merge(b_L_sub_sub,h_L_sub_sub,by="SUBJID",all.x=TRUE)
b_L_sub_sub1$Study <-"LIGHT"

dim(b_L_sub_sub1)

b_all_sub_sub1<-rbind(b_all_sub_sub,b_L_sub_sub1)

####Need a separate file to get information on histology, location and FIGO stage

###



####now merge with biomarker data
colnames(b_all_sub_sub1)<-c("SUBJECT","Race","Histology","Location","FIGO Stage","Study")
b_all_sub_sub1$key<-paste(b_all_sub_sub1$SUBJECT,b_all_sub_sub1$Study,sep="_")
all_HRD$key<-paste(all_HRD$SUBJECT,all_HRD$Study,sep="_")

link<-intersect(all_HRD$key,b_all_sub_sub1$key)

all_HRD_base<-merge(all_HRD,b_all_sub_sub1,by="key",all=TRUE)
write.csv(all_HRD_base,file="/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_2024_0529_with_base_char.csv",row.names=FALSE)

all_HRD_base<-read.csv("/Documents/Clinical_studies/Lynparza/Olaparib_Ovarian/Cross_study/all_HRD_2024_0529_with_base_char.csv")
all_HRD_base_sub<-all_HRD_base[!is.na(all_HRD_base$Study.x),]
all_HRD_base_sub<-all_HRD_base_sub[!is.na(all_HRD_base_sub$Study.y),]
all_HRD_base_sub<-all_HRD_base_sub[!all_HRD_base_sub$HRD_status =="HRD unknown (not tested/test failed)",]
dim(all_HRD_base_sub)
table(all_HRD_base_sub$Study.x)

####RACE
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="PAOLA"])
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="OPINION"])
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="LIGHT"])
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="Study 19"])
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="SOLO1" ])
table(all_HRD_base_sub$Race[all_HRD_base_sub$Study.x =="SOLO2" ])

table(all_HRD_base_sub$Race[all_HRD_base_sub$HRD_status !="HRD unknown (not tested/test failed)"])

cut<-c("FAILED",NA)
summary(as.numeric(as.character(all_HRD_base_sub$HRD_Score[!all_HRD_base_sub$HRD_Score %in% cut]))) ###This includes PAOLA-1
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00   31.00   54.00   49.16   66.00   96.00 

all_HRD_base_sub_GIS<-all_HRD_base_sub[!all_HRD_base_sub$HRD_Score %in% cut,]

table(all_HRD_base_sub_GIS$Race)


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="White"])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00   38.00   57.00   52.01   67.00   96.00    

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="White" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Black or African American"])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00   41.00   53.00   49.50   62.25   86.00 

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Black or African American" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.00   52.00   59.00   57.79   68.25   92.00



summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Asian"])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.00   52.00   59.00   57.79   68.25   92.00 
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Asian" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Other/Unknown"])))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#18.00   38.00   46.50   46.19   53.50   70.00 
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Race =="Other/Unknown" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[!all_HRD_base_sub_GIS$Race =="NA"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[!all_HRD_base_sub_GIS$Race =="NA" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))



table(all_HRD_base_sub_GIS$Status_3[!all_HRD_base_sub_GIS$Race =="NA"])

#Non-tBRCAm     tBRCAm 
#422        541 

###56.2% (541/963) of GIS evaluable population are BRCAm


table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Race =="White"])

#Non-tBRCAm     tBRCAm 
#384        463 

#54.6% (463/847) of GIS evaluable population that are White are BRCAm

table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Race =="Black or African American"])
#Non-tBRCAm     tBRCAm 
#11         13 

#54.2% (13/24) of GIS evaluable population that are White are BRCAm

table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Race =="Asian"])
#Non-tBRCAm     tBRCAm 
#21         55 
#72.4% (55/76) of GIS evaluable population that are White are BRCAm

table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Race =="Other/Unknown"])
#Non-tBRCAm     tBRCAm 
#6         10 
#62.5% (10/16) of GIS evaluable population that are White are BRCAm


###Now make the plot for Race

all_HRD_base_sub_GIS_noNA<-all_HRD_base_sub_GIS[!is.na(all_HRD_base_sub_GIS$Race),]
data.frame(table(all_HRD_base_sub_GIS_noNA$Race))
data.frame(table(all_HRD_base_sub_GIS_noNA$Race[all_HRD_base_sub_GIS_noNA$Status_3 =="tBRCAm"]))

all_HRD_base_sub_GIS_noNA$Race<-factor(all_HRD_base_sub_GIS_noNA$Race, levels=c("Other/Unknown","Black or African American","Asian","White"))
all_HRD_base_sub_GIS_noNA$HRD_Score <-as.numeric(as.character(all_HRD_base_sub_GIS_noNA$HRD_Score))

ggplot(all_HRD_base_sub_GIS_noNA, aes(x = Race , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Race),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=18))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


t.test(all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="White"],all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Black or African American"]) ###NS
t.test(all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="White"],all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Asian"]) ###p-value = 0.007968 (higher proportion of BRCAm is Asian pop)
t.test(all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="White"],all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Other/Unknown"]) ###NS
t.test(all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Asian"],all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Black or African American"]) ###NA
t.test(all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Asian"],all_HRD_base_sub_GIS_noNA$HRD_Score[all_HRD_base_sub_GIS_noNA$Race =="Other/Unknown"]) ###p-value = 0.01176 (higher proportion of BRCAm in Asian pop)


all_HRD_base_sub_GIS_noNA_BRCA<-subset(all_HRD_base_sub_GIS_noNA,all_HRD_base_sub_GIS_noNA$Status_3 =="tBRCAm")

ggplot(all_HRD_base_sub_GIS_noNA_BRCA, aes(x = Race , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Race),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))


t.test(all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="White"],all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Black or African American"])
t.test(all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="White"],all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Asian"]) ##NS
t.test(all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="White"],all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Other/Unknown"])
t.test(all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Asian"],all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Black or African American"])
t.test(all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Asian"],all_HRD_base_sub_GIS_noNA_BRCA$HRD_Score[all_HRD_base_sub_GIS_noNA_BRCA$Race =="Other/Unknown"]) 









####HISTOLOGY
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="PAOLA"])
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="OPINION"])
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="LIGHT"])
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="Study 19"])
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="SOLO1"])
table(all_HRD_base_sub$Histology[all_HRD_base_sub$Study.x =="SOLO2"])
table(all_HRD_base_sub$Histology)
table(all_HRD_base_sub_GIS$Histology)
table(all_HRD_base_sub_GIS$Status_3)

table(all_HRD_base_sub_GIS$Histology[all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score)))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Serous"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Serous" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Endometrioid"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Endometrioid" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Other/Unknown"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Other/Unknown" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

table(all_HRD_base_sub_GIS$Status_3)
table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Histology =="Serous"])
table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Histology =="Endometrioid"])
table(all_HRD_base_sub_GIS$Status_3[all_HRD_base_sub_GIS$Histology =="Other/Unknown"])


#### Now make plot for Histology
all_HRD_base_sub_GIS$HRD_Score <-as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score))
all_HRD_base_sub_GIS$Histology<-factor(all_HRD_base_sub_GIS$Histology, levels=c("Other/Unknown","Endometrioid", "Serous"))
table(all_HRD_base_sub_GIS$Histology)

ggplot(all_HRD_base_sub_GIS, aes(x = Histology , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Histology),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Serous"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Endometrioid"]) ###NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Serous"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Other/Unknown"]) ###NS)
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Endometrioid"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Histology =="Other/Unknown"]) ###NS



all_HRD_base_sub_GIS_BRCA<-subset(all_HRD_base_sub_GIS,all_HRD_base_sub_GIS$Status_3 =="tBRCAm")
table(all_HRD_base_sub_GIS_BRCA$Histology)

ggplot(all_HRD_base_sub_GIS_BRCA, aes(x = Histology , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Histology),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))



t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Serous"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Endometrioid"]) ###NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Serous"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Other/Unknown"]) ###NS)
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Endometrioid"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Histology =="Other/Unknown"]) ###NS





###FIGO
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="PAOLA" ])
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="OPINION" ])
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="LIGHT" ])
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="Study 19"])
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="SOLO1"])
table(all_HRD_base_sub$FIGO.Stage[all_HRD_base_sub$Study.x =="SOLO2"])
table(all_HRD_base_sub$FIGO.Stage)
table(all_HRD_base_sub_GIS$FIGO.Stage)
table(all_HRD_base_sub_GIS$FIGO.Stage[all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="II"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="II" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="III"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="III" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="IV"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="IV" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))


summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="Other/Unknown"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="Other/Unknown" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))




#### Now make plot for FIGO
all_HRD_base_sub_GIS$FIGO.Stage<-factor(all_HRD_base_sub_GIS$FIGO.Stage, levels=c("Other/Unknown","IV","III","II"))


ggplot(all_HRD_base_sub_GIS, aes(x = FIGO.Stage , y = HRD_Score),color="black")+
  geom_violin(aes(fill=FIGO.Stage),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="II"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="III"]) ###NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="II"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="IV"]) ###NS)
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="II"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="Other/Unknown"]) ###NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="III"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="Other/Unknown"]) ###p-value = NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="IV"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$FIGO.Stage =="Other/Unknown"]) ###NS



all_HRD_base_sub_GIS_BRCA<-subset(all_HRD_base_sub_GIS,all_HRD_base_sub_GIS$Status_3 =="tBRCAm")

ggplot(all_HRD_base_sub_GIS_BRCA, aes(x = FIGO.Stage , y = HRD_Score),color="black")+
  geom_violin(aes(fill=FIGO.Stage),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))



t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="II"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="III"]) ###NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="II"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="IV"]) ###NS)
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="II"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="Other/Unknown"]) ###NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="III"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="Other/Unknown"]) ###NS when consider BRCAm
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="IV"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$FIGO.Stage =="Other/Unknown"]) ###NS when consider BRCAm




###LOCATION
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="PAOLA" ])
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="OPINION"])
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="LIGHT"])
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="Study 19"])
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="SOLO1"])
table(all_HRD_base_sub$Location[all_HRD_base_sub$Study.x =="SOLO2"])
table(all_HRD_base_sub$Location)
table(all_HRD_base_sub_GIS$Location)
table(all_HRD_base_sub_GIS_BRCA$Location)

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score)))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Ovary"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Ovary" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Fallopian Tube"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Fallopian Tube" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Peritoneal"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Peritoneal" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Other/Unknown"])))
summary(as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Other/Unknown" & all_HRD_base_sub_GIS$Status_3 =="tBRCAm"])))

#### Now make plot for Location
#all_HRD_base_sub_GIS$HRD_Score <-as.numeric(as.character(all_HRD_base_sub_GIS$HRD_Score))
#all_HRD_base_sub_GIS$Location<-factor(all_HRD_base_sub_GIS$Location, levels=c("Other/Unknown","Fallopian Tube or Primary Peritoneal","Ovary"))
all_HRD_base_sub_GIS$Location<-factor(all_HRD_base_sub_GIS$Location, levels=c("Other/Unknown","Peritoneal","Fallopian Tube","Ovary"))


ggplot(all_HRD_base_sub_GIS, aes(x = Location , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Location),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))

t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Ovary"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Fallopian Tube"]) ### NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Ovary"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Peritoneal"]) ### 
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Ovary"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Other/Unknown"]) ###NS)
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Fallopian Tube"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Peritoneal"]) ###NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Peritoneal"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Other/Unknown"]) ###NS
t.test(all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Fallopian Tube"],all_HRD_base_sub_GIS$HRD_Score[all_HRD_base_sub_GIS$Location =="Other/Unknown"]) ###NS



all_HRD_base_sub_GIS_BRCA<-subset(all_HRD_base_sub_GIS,all_HRD_base_sub_GIS$Status_3 =="tBRCAm")

ggplot(all_HRD_base_sub_GIS_BRCA, aes(x = Location , y = HRD_Score),color="black")+
  geom_violin(aes(fill=Location),size=1)+
  xlab("")+ 
  coord_flip()+
  geom_jitter(color="black",shape=16, position=position_jitter(0.2))+
  geom_hline(yintercept=42,linetype="dashed")+
  #facet_grid(.~ Status_3,scales="free")+
  scale_y_continuous(limits = c(-1, 101), breaks = seq(0, 101, by = 25))+
  #scale_fill_manual(values=c("gBRCA1m"="blueviolet","gBRCA2m"="blueviolet","sBRCA1m"="cyan3","sBRCA2m"="cyan3",`s/g BRCA not determined BRCA1m`= "grey50",`s/g BRCA not determined BRCA2m`= "grey50",`HRRm and non-BRCAm` ="orange","Non-HRRm" ="forestgreen"))+
  theme_bw()+
  theme(strip.text.x = element_text(size=18, color="Black",
                                    face="bold.italic"))+
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid"))+
  theme(axis.text.x=element_text(size=25),axis.title.x=element_text(size=25,face="bold"))+
  theme(axis.text.y=element_text(size=25),axis.title.y =element_text(size=25,face="bold"))+
  theme(legend.title=element_text(size=19,face="bold"))+
  theme(legend.text=element_text(size=20))+
  ylab("HRD Score")+
  theme(legend.position = 'top')+
  guides(fill=guide_legend(title=""))



t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Ovary"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Fallopian Tube"]) ### NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Ovary"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Peritoneal"]) ### NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Ovary"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Other/Unknown"]) ###NS)
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Fallopian Tube"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Peritoneal"]) ###NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Peritoneal"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Other/Unknown"]) ###NS
t.test(all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Fallopian Tube"],all_HRD_base_sub_GIS_BRCA$HRD_Score[all_HRD_base_sub_GIS_BRCA$Location =="Other/Unknown"]) ###NS






