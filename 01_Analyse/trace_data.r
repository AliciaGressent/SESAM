# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/01_Analyse/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(raster)
library(maptools)
library(RColorBrewer)
library(fields)
library(rgdal)
library(ggplot2)
library(chron)

#################################################################################
#                               DATA VISUALIZATION                              #
#                                Read sensor data                               #
#       plot time series and calibration plot at reference stations             #
#                Created 05/06/2018 and updated 13/05/2020                      #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#             FUNCTIONS             #
#####################################

# Function to average sensor repilicas at reference station
avgsens <- function(i,data1_sta_all,data2_sta_all){
    data1_sta <- subset(stagrid, stagrid$id_sensor == nsensor[i]) # Select sensor data
    data1_sta_all <- rbind(data1_sta_all,data1_sta)
    data1_sta$interval_time <- cut(data1_sta$datetime,breaks=timevec_out)
    data2_sta <- aggregate(. ~ interval_time,data=data1_sta,mean,na.action=na.pass) # Aggregate data over a 15 minutes period
    data2_sta$datetime = as.POSIXct(as.character(data2_sta$interval_time),origin = "1970-01-01",tz='GMT')
    data2_sta_all = rbind(data2_sta_all,data2_sta) #concatenate for nsensor at the station
    newlist <- list(data1_sta_all,data2_sta_all)
    return(newlist)
}

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/INPUTS/" # path for input directory
dirout <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/01_Analyse/figs/" # path for output directory
start_date <-"2018-11-01 00:00:00" # starting date
end_date <-"2018-11-30 23:59:59" # ending date
pol <- "PM10" # pollutant
ref_file <- paste0("mesures_ref_qthourly_",pol,".txt") # reference measurements file name
Fsens_file <- "dataout_fixe_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
sta1_name <- "Bouteillerie" # 1st reference station's name
sta2_name <- "Victor Hugo" # 2nd reference station's name
code_sta_pol <- c(188,140,238,107,239) # reference station's code
lat1 = 47.223253; lon1 = -1.539075 # Station N°1
lat2 = 47.204264; lon2 = -1.552927 # Station N°2
lat3 = 47.194049; lon3 = -1.581454 # Station N°3
lat4 = 47.252721; lon4 = -1.574653 # Station N°4
lat5 = 47.185581; lon5 = -1.591332 # Station N°5

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
nsta_pol <- unique(refpol$code_station) # Number of stations
lat_refpol <- rep(0,length(refpol[,1])); lon_refpol <- rep(0,length(refpol[,1]))
lon_refpol[which(refpol$code_station==code_sta_pol[1])]=lon1; lat_refpol[which(refpol$code_station==code_sta_pol[1])]=lat1;
lon_refpol[which(refpol$code_station==code_sta_pol[2])]=lon2; lat_refpol[which(refpol$code_station==code_sta_pol[2])]=lat2;
lon_refpol[which(refpol$code_station==code_sta_pol[3])]=lon3; lat_refpol[which(refpol$code_station==code_sta_pol[3])]=lat3;
lon_refpol[which(refpol$code_station==code_sta_pol[4])]=lon4; lat_refpol[which(refpol$code_station==code_sta_pol[4])]=lat4;
lon_refpol[which(refpol$code_station==code_sta_pol[5])]=lon5; lat_refpol[which(refpol$code_station==code_sta_pol[5])]=lat5;
refpol <- cbind(refpol,lat_refpol,lon_refpol)

#####################################
#         READ SENSOR DATA          #
#####################################

# Read fixed sensor mearsurements of the pollutant
staFile <- paste0(indir,Fsens_file)
stagrid <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
stagrid$datetime = as.POSIXct(as.character(stagrid$datetime),tz='GMT')
colnames(stagrid)[which(names(stagrid) == tolower(pol))] <- "pol"
nsensor <- unique(stagrid$id_sensor) # Number of sensors
sta <- refpol[which(refpol$code_station==code_sta_pol[2]),] # Create a time vector with intervals defined by the resolution of the reference station
timevec_out <- sta$date

####################################
#    AVG SENSOR DATA AT STATIONS   # 
####################################

xmin=as.POSIXct(start_date,origin = "1970-01-01",tz = "GMT") # Set time limits
xmax=as.POSIXct(end_date,origin = "1970-01-01",tz = "GMT")

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

data3_sta1_all <- aggregate(.~lat+lon+datetime,data=data2_sta1_all,mean,na.action=na.omit) # AVG of the 3 sensors
data3_sta2_all <- aggregate(.~lat+lon+datetime,data=data2_sta2_all,mean,na.action=na.omit) # AVG of the 3 sensors

# When observations are missing, select the date for which data exist
dataout_sta1 <- c()
for (l in 1:length(ref_sta1_pol_sub[,1])){
        tmp = data3_sta1_all[which(data3_sta1_all$datetime==ref_sta1_pol_sub$date[l]),]
        if (length(tmp[,1] ) == 0){tmp = data3_sta1_all[1,]+NA }
        dataout_sta1 <- rbind(dataout_sta1, tmp)
}
dataout_sta2 <- data3_sta2_all

# Derive correlation 
dataout_sta1$pol_ref <- ref_sta1_pol_sub$pol
dataout_sta1 = na.omit(dataout_sta1)
correlation_sta1 = round(cor(dataout_sta1$pol_ref,dataout_sta1$pol),2)

dataout_sta2$pol_ref <- ref_sta2_pol_sub$pol
dataout_sta2 = na.omit(dataout_sta2)
correlation_sta2 = round(cor(dataout_sta2$pol_ref,dataout_sta2$pol),2)

####################################
#               PLOT               #
####################################

# Station N°1
df1 <- subset(data3_sta1_all,select=-c(lon,lat,interval_time,id_sensor,temperature,humidity,pressure,pm25))
var1 <- as.vector(t(rep("Fixed LCS avg 15 min",length(data3_sta1_all[,1]))))
df1$VAR <- var1

df2 <- subset(ref_sta1_pol_sub,select=-c(lat_refpol,lon_refpol,code_station))
var2 <- as.vector(t(rep("Reference",length(ref_sta1_pol_sub[,1]))))
df2$VAR <- var2
names(df2)[names(df2) == "date"] <- "datetime"

df3 <- rbind(df1,df2)
df3$VAR = as.factor(df3$VAR)

# Station N°2
df4 <- subset(data3_sta2_all,select=-c(lon,lat,interval_time,id_sensor,temperature,humidity,pressure,pm25))
var4 <- as.vector(t(rep("Fixed LCS avg 15 minLCS",length(data3_sta2_all[,1]))))
df4$VAR <- var4

df5 <- subset(ref_sta2_pol_sub,select=-c(lat_refpol,lon_refpol,code_station))
var5 <- as.vector(t(rep("Reference",length(ref_sta2_pol_sub[,1]))))
df5$VAR <- var5
names(df5)[names(df5) == "date"] <- "datetime"

df6 <- rbind(df4,df5)
df6$VAR = as.factor(df6$VAR)

# Do plot
source('multiplot.r')
png(filename=paste0(dirout,"Figure1.png"), width=1200, height=800, type="cairo",bg = "white")
p1 <- ggplot(df6, aes(x=datetime, y=pol, colour=VAR)) +
geom_line(size=1,alpha=0.7)+ylim(0,200)+
scale_color_manual(values=c("#006699","gray57"))+
labs(title=sta2_name,y=bquote(.(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.title.x = element_blank(),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p2 <- ggplot(df3, aes(x=datetime, y=pol, colour=VAR)) +
geom_line(size=1,alpha=0.7)+ylim(0,200)+
scale_color_manual(values=c("#006699","gray57"))+
labs(title=sta1_name,y=bquote(.(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=16),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        axis.title.x = element_blank(),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p3 <- ggplot(dataout_sta1, aes(x=pol_ref, y=pol)) +
geom_point(size=1,color="black")+
geom_smooth(method=lm, se=FALSE,color="red")+
geom_abline(intercept = 0,linetype="dashed", color = "red",size=1)+
annotate("text", x=40 , y=130, label = paste0("R=",correlation_sta1),size=6)+
ylim(0,140)+xlim(0,220)+
labs(title="",y=bquote(LCS ~ .(pol) ~ (mu*g/m^3)),x=bquote(Reference ~ .(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p4 <- ggplot(dataout_sta2, aes(x=pol_ref, y=pol)) +
geom_point(size=1,color="black")+
geom_smooth(method=lm, se=FALSE, color="red")+
geom_abline(intercept = 0,linetype="dashed", color = "red",size=1)+
annotate("text", x=40 , y=130, label = paste0("R=",correlation_sta2),size=6)+
ylim(0,140)+xlim(0,220)+
labs(title="",y=bquote(LCS ~ .(pol) ~ (mu*g/m^3)),x=bquote(Reference ~ .(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        legend.text = element_text(size =14),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
pf <- multiplot(p1,p3,p2,p4,layout=matrix(c(1,1,4,3,3,2), 2, 3, byrow = TRUE),cols=2,fontsize = 24)
pf
dev.off()


