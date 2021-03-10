# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/01_Stat/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(raster)
library(maptools)
library(RColorBrewer)
library(fields)
library(rgdal)
library(ggplot2)
library(chron)

#################################################################################
#                             DATA VISUALIZATION                                #
#                   Read sensor data and plot histogams                         #
#                 Created 22/10/2019 and updated 13/05/2020                     #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#             FUNCTIONS             #
#####################################

# Function to average sensor repilicas at reference station
avgsens <- function(i,data1_sta_all,data2_sta_all){
    data1_sta <- subset(stagrid_fs, stagrid_fs$id_sensor == nsensor[i]) # Select sensor data
    data1_sta_all <- rbind(data1_sta_all,data1_sta)
    data1_sta$interval_time <- cut(data1_sta$datetime,breaks=timevec_out)
    data2_sta <- aggregate(. ~ interval_time,data=data1_sta,mean,na.action=na.pass) # Aggregate data over a 15 minutes period
    data2_sta$datetime = as.POSIXct(as.character(data2_sta$interval_time),origin = "1970-01-01",tz='GMT')
    data2_sta_all = rbind(data2_sta_all,data2_sta) # Concatenate for nsensor at the station
    newlist <- list(data1_sta_all,data2_sta_all)
    return(newlist)
}

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/INPUTS/" # path for input directory
dirout <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/01_Stat/figs/" # path for output directory
pol="PM10" # pollutant
ref_file <- paste0("mesures_ref_qthourly_",pol,".txt") # reference measurements file name
Fsens_file <- "dataout_fixe_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
Msens_file <- "dataout_mobile_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
code_sta_pol <- c(188,140,238,107,239) # reference station's code
sta_name <- c("Bouteillerie","Victor Hugo","Trentemoult","Chavinière","Les Coûets") # reference station's name

#####################################
#           READ REF DATA           #
#####################################

# Read reference station measurements of the pollutant
refpol <- read.table(file=paste0(indir,ref_file), header=TRUE, sep="\t",skip=0)
refpol <- subset(refpol, select=-c(nom_station,polluant))
colnames(refpol)[which(names(refpol) == "niveau")] <- "pol"
refpol$date <- gsub("T", " ", refpol$date)
refpol$date = as.POSIXct(as.character(refpol$date),tz='GMT')
refpol=refpol[!duplicated(refpol),] # Remove duplicated data
refpol <- na.omit(refpol) # Remove NAN
nsta_pol <- unique(refpol$code_station)

# Create a time vector with intervals defined by the resolution of the ref station
sta <- refpol[which(refpol$code_station==code_sta_pol[1]),]
timevec_out <- sta$date

#####################################
#         READ SENSOR DATA          #
#####################################

staFile <- paste0(indir,Fsens_file)
stagrid_fs <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
stagrid_fs$datetime = as.POSIXct(as.character(stagrid_fs$datetime),tz='GMT')
colnames(stagrid_fs)[which(names(stagrid_fs) == tolower(pol))] <- "pol"
nsensor <- unique(stagrid_fs$id_sensor) # Sensors ID
n1 <- length(nsensor) # Number of sensors
n2 <- n1/2 # Half the number of sensor

staFile <- paste0(indir,Msens_file)
stagrid_ms <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
stagrid_ms$datetime = as.POSIXct(as.character(stagrid_ms$datetime),tz='GMT')
names(stagrid_ms)[9]<- 'pol'
nsensorb <- unique(stagrid_ms$id_sensor) # Sensors ID
n1b <- length(nsensorb) # Number of sensors
n2b <- ceiling(n1b/2) # Half the number of sensor

#####################################
#     COLOCATED FS AND REF DATA     #
#####################################

xmin=as.POSIXct("2018-11-01 00:00:00",origin = "1970-01-01",tz = "GMT") # Set time limits
xmax=as.POSIXct("2018-11-30 23:59:59",origin = "1970-01-01",tz = "GMT")

# Init dataframes
data1_sta_all <- c(); data2_sta_all <- c()

for (i in 1:3){ tmp=avgsens(i,data1_sta_all,data2_sta_all)
data1_sta_all=as.data.frame(tmp[1]); data2_sta_all=as.data.frame(tmp[2])
}
data2_sta1_all <- data2_sta_all

data1_sta_all <- c(); data2_sta_all <- c()
for (i in 4:6){ tmp=avgsens(i,data1_sta_all,data2_sta_all)
data1_sta_all=as.data.frame(tmp[1]);data2_sta_all=as.data.frame(tmp[2])
}
data2_sta2_all <- data2_sta_all

ref_sta1_pol <- subset(refpol,refpol$code_station == code_sta_pol[1]) # Select reference measurements station N°1
ref_sta1_pol_sub <-  subset(ref_sta1_pol, ref_sta1_pol$date >= xmin & ref_sta1_pol$date <= xmax) # Subset data ref for the estimation period
ref_sta2_pol <- subset(refpol,refpol$code_station == code_sta_pol[2]) # Select reference measurements station N°2
ref_sta2_pol_sub <-  subset(ref_sta2_pol, ref_sta2_pol$date >= xmin & ref_sta2_pol$date <= xmax) # Subset data ref for the estimation period

data_fs_sta1_all <- aggregate(.~lat+lon+datetime,data=data2_sta1_all,mean,na.action=na.omit) # AVG of the 3 sensors
data_fs_sta2_all <- aggregate(.~lat+lon+datetime,data=data2_sta2_all,mean,na.action=na.omit) # AVG of the 3 sensors

#####################################
#       PREPARE DATA FOR PLOT       #
#####################################

# Subset data for the time period
refpol_sub <- subset(refpol, refpol$date >= xmin & refpol$date <= xmax)
stagrid_fs_sub <- subset(stagrid_fs, stagrid_fs$datetime >= xmin & stagrid_fs$datetime <= xmax)
stagrid_ms_sub <- subset(stagrid_ms, stagrid_ms$datetime >= xmin & stagrid_ms$datetime <= xmax)

# Reference and FS avg 15 min data colocated at station N°2
pol_ref = as.numeric(ref_sta2_pol_sub$pol)
pol_fs_avg = data_fs_sta2_all$pol
var1 <- as.vector(t(rep("Reference",length(pol_ref))))
var2 <- as.vector(t(rep("Fixed LCS avg 15min",length(pol_fs_avg))))
df1 <- as.data.frame(var1)
df1$pol <- pol_ref
names(df1)[names(df1) == "var1"] <- "VAR"
df2 <- as.data.frame(var2)
df2$pol <- pol_fs_avg
names(df2)[names(df2) == "var2"] <- "VAR"
df <- rbind(df2,df1)

# Reference and FS avg 15 min data colocated at station N°1
pol_refb = as.numeric(ref_sta1_pol_sub$pol)
pol_fsb_avg = data_fs_sta1_all$pol
var1 <- as.vector(t(rep("Reference",length(pol_refb))))
var2 <- as.vector(t(rep("Fixed LCS avg 15min",length(pol_fsb_avg))))
df1b <- as.data.frame(var1)
df1b$pol <- pol_refb
names(df1b)[names(df1b) == "var1"] <- "VAR"
df2b <- as.data.frame(var2)
df2b$pol <- pol_fsb_avg
names(df2b)[names(df2b) == "var2"] <- "VAR"
dfb <- rbind(df2b,df1b)

# FS and MS data no avg
pol_fs = stagrid_fs_sub$pol
pol_ms = stagrid_ms_sub$pol
pol_ref_all <- as.numeric(refpol_sub$pol)
var3 <- as.vector(t(rep("Fixed LCS",length(pol_fs))))
var4 <- as.vector(t(rep("Mobile LCS",length(pol_ms))))
var5 <- as.vector(t(rep("Reference",length(pol_ref_all))))
df3 <- as.data.frame(var3)
df3$pol <- pol_fs
names(df3)[names(df3) == "var3"] <- "VAR"
df4 <- as.data.frame(var4)
df4$pol <- pol_ms
names(df4)[names(df4) == "var4"] <- "VAR"
df5 <- as.data.frame(var5)
df5$pol <- pol_ref_all
names(df5)[names(df5) == "var5"] <- "VAR"
dfc <- rbind(df3,df4,df5)

# Datetime ref, FS and FS
date_ref = refpol_sub$date
date_fs = stagrid_fs_sub$datetime
date_ms = stagrid_ms_sub$datetime

var1 <- as.vector(t(rep("Reference",length(date_ref))))
df1_date = as.data.frame(var1)
df1_date$date <- date_ref
names(df1_date)[names(df1_date) == "var1"] <- "VAR"

var3 <- as.vector(t(rep("Fixed LCS",length(date_fs))))
df3_date = as.data.frame(var3)
df3_date$date <- date_fs
names(df3_date)[names(df3_date) == "var3"] <- "VAR"

var4 <- as.vector(t(rep("Mobile LCS",length(date_ms))))
df4_date = as.data.frame(var4)
df4_date$date <- date_ms
names(df4_date)[names(df4_date) == "var4"] <- "VAR"

df_date <- rbind(df3_date,df4_date,df1_date)

#####################################
#            DO PLOT                #
#####################################

source('multiplot.r')
png(filename=paste0(dirout,"Figure2.png"), width=800, height=800, type="cairo",bg = "white")

p1 <- ggplot(df, aes(x=pol, fill=VAR, color=VAR)) +
geom_density(alpha=0.5)+
xlim(c(0, 200))+ylim(c(0, 0.05))+
labs(title="a)",x=bquote(.(pol) ~ (mu*g/m^3)), y = "Density")+
annotate("text", x=150 , y=0.045, label = "Victor Hugo",size=6)+
scale_color_manual(values=c("#006699","#666666"))+
scale_fill_manual(values=c("#006699","#666666"))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p1
p2 <- ggplot(dfb, aes(x=pol, fill=VAR, color=VAR)) +
geom_density(alpha=0.5)+
xlim(c(0, 200))+ ylim(c(0, 0.05))+
labs(title="b)",x=bquote(.(pol) ~ (mu*g/m^3)), y = "Density")+
annotate("text", x=150 , y=0.045, label = "Bouteillerie",size=6)+
scale_color_manual(values=c("#006699","#666666"))+
scale_fill_manual(values=c("#006699","#666666"))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p2
p3 <- ggplot(dfc, aes(x=pol, fill=VAR, color=VAR)) +
geom_density(alpha=0.5)+
xlim(c(0, 200))+ylim(c(0, 0.05))+
labs(title="c)",x=bquote(.(pol) ~ (mu*g/m^3)), y = "Density")+
scale_color_manual(values=c("#006699","#CC6600","#666666"))+
scale_fill_manual(values=c("#006699","#CC6600","#666666"))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p3
p4 <- ggplot(df_date, aes(x=date, fill=VAR, color=VAR)) +
geom_density(alpha=0.5)+
labs(title="d)",x="Sampling time", y = "Density")+
scale_color_manual(values=c("#006699","#CC6600","#666666"))+
scale_fill_manual(values=c("#006699","#CC6600","#666666"))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p4
pf <- multiplot(p1,p3,p2,p4,cols=2,fontsize = 24)
pf
dev.off()

