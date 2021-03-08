# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/01_Stat/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(raster)
library(maptools)
library(RColorBrewer)
library(fields)
library(rgdal)
library(ggplot2)
library(chron)

#####################################################################################
#           CALCULATE REPEATABILITY AT LOW CONCENTRATION FOR SENSORS                #
#                           From Spinelle et al., 2013                              #
#     Identify a period with no PM variation in reference data and compute Rep      #
#                Created 05/04/2019 and updated on 14/05/2020                       #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr             #
#####################################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/INPUTS/" # path for input directory
dirout <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/01_Stat/figs/" # path for output directory
pol="PM10"
Fsens_file <- "dataout_fixe_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
ref_file <- paste0("mesures_ref_qthourly_",pol,".txt") # reference measurements file name
code_sta_pol <- c(188,140,238,107,239) # reference station's code
sta_name <- c("Bouteillerie","Victor Hugo") # reference station's name
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
#          SELECT PERIODS           #
#####################################

# Init dataframes
len_dataout <- c()
dataout1 <- c();dataout2 <- c();dataout3 <- c();dataout4 <- c()
dataout5 <- c();dataout6 <- c();dataout7 <- c();dataout8 <- c()
dataout9 <- c();dataout10 <- c();dataout11 <- c();dataout12 <- c()
dataout13 <- c();dataout14 <- c();dataout15 <- c();dataout16 <- c()
P1 <- c(); P2 <- c(); P3 <- c()

# Select data at station
pp=2 # Choose the station from which you want to select the period, station N°2 !!! MUST BE UPDATED BY THE USER !!!
refpol_sta <- refpol[which(refpol$code_station==code_sta_pol[2]),] # !!! MUST BE UPDATED BY THE USER !!! 

# Look for period of at least 30 minutes for which there is no PM variation for low pollutant concentrations ([PM10]=3ug/m3)
refpol_sta_low=refpol_sta[which(refpol_sta$pol==3),]

# Init variables for the loop
dataout <- c(); n=0; m=0

if (pol == "PM10"){
    output_var <- c("dataout1","dataout2","dataout3","dataout4","dataout5","dataout6","dataout7","dataout8","dataout9","dataout10","dataout11","dataout12","dataout13","dataout14","dataout15","dataout16")
}else if (pol== "PM25"){
    output_var <- c("dataout1","dataout2","dataout3","dataout4","dataout5","dataout6","dataout7","dataout8","dataout9","dataout10","dataout11","dataout12","dataout13","dataout14","dataout15","dataout16",
                    "dataout17","dataout18","dataout19","dataout20","dataout21","dataout22","dataout23","dataout24","dataout25","dataout26","dataout27","dataout28","dataout29","dataout30","dataout31","dataout32",
                    "dataout33","dataout34","dataout35","dataout36","dataout37","dataout38","dataout39","dataout40","dataout41","dataout42","dataout43","dataout44","dataout45","dataout46","dataout47","dataout48")
}

new_period="TRUE"

# Separate data period
for (i in 1:100){
    m=m+1
    if (new_period=="TRUE" & n>0){m=m-1}
    if (new_period=="TRUE"){
        dataout<-refpol_sta_low[m,]
        new_period="FALSE"
    }else{
        date1 = refpol_sta_low$date[m-1]
        date2 = refpol_sta_low$date[m]
        diff_date = difftime(date2,date1,units="hours")
        diff_date = as.numeric(diff_date)
        if (diff_date == 0.25){
            dataout <- rbind(dataout,refpol_sta_low[m,])
        }else{
            new_period="TRUE"
            n=n+1
            assign(output_var[n], dataout)
        }
    }
    if (m==length(refpol_sta_low[,1])){break}
}

# Create a vector of length of the dataout periods

if (pol == "PM10"){
    len_dataout <- c(length(dataout1[,1]),length(dataout2[,1]),length(dataout3[,1]),length(dataout4[,1]),length(dataout5[,1]),length(dataout6[,1]),
                    length(dataout7[,1]),length(dataout8[,1]),length(dataout9[,1]),length(dataout10[,1]),length(dataout11[,1]),length(dataout12[,1]))
}else if (pol== "PM25"){
    len_dataout <- c(length(dataout1[,1]),length(dataout2[,1]),length(dataout3[,1]),length(dataout4[,1]),length(dataout5[,1]),length(dataout6[,1]),
                    length(dataout7[,1]),length(dataout8[,1]),length(dataout9[,1]),length(dataout10[,1]),length(dataout11[,1]),length(dataout12[,1]),
                    length(dataout13[,1]),length(dataout14[,1]),length(dataout15[,1]),length(dataout16[,1]),length(dataout17[,1]),length(dataout18[,1]),
                    length(dataout19[,1]),length(dataout20[,1]),length(dataout21[,1]),length(dataout22[,1]),length(dataout23[,1]),length(dataout24[,1]),
                    length(dataout25[,1]),length(dataout26[,1]),length(dataout27[,1]))
}

idx=which(len_dataout>=5) # Select periods with > 5 points

# Assign period name to each subset data
ref_period <- c("P1","P2","P3")
for (n in 1:length(idx)){
    tmp=eval(parse(text=paste0("dataout",idx[n])))
    tmp2=tmp[1:5,]
    assign(ref_period[n], tmp2)
}

#Subset for P1
tdebP1 <- min(P1$date); tendP1 <- max(P1$date)
tdebP2 <- min(P2$date); tendP2 <- max(P2$date)
tdebP3 <- min(P3$date); tendP3 <- max(P3$date)

#####################################
#        SELECT SENSOR DATA         #
#####################################

# Read fixed sensor mearsurements
staFile <- paste0(indir,Fsens_file)
stagrid <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
stagrid$datetime = as.POSIXct(as.character(stagrid$datetime),tz='GMT')
colnames(stagrid)[which(names(stagrid) == tolower(pol))] <- "pol"
nsensor <- unique(stagrid$id_sensor) # Number of sensors

# Init vector
Rep_mean_all_sensors <- c()

for (i in 4:6) { # ith sensors corresponding to the station !!! MUST BE UPDATED BY THE USER !!!

  # Select sensor data and average depending on the deltatime
  FS <- subset(stagrid, stagrid$id_sensor == nsensor[i])
  FS_P1 <- subset(FS, FS$datetime >= tdebP1 & FS$datetime <= tendP1)
  FS_P2 <- subset(FS, FS$datetime >= tdebP2 & FS$datetime <= tendP2)
  FS_P3 <- subset(FS, FS$datetime >= tdebP3 & FS$datetime <= tendP3)
    
  FS_P1_moy=mean(FS_P1$pol)
  Sr_P1 = sqrt((sum((FS_P1$pol - FS_P1_moy)**2))/(length(FS_P1[,1])-1))    
  Rep_P1 = 2*sqrt(2*Sr_P1)

  FS_P2_moy=mean(FS_P2$pol)
  Sr_P2 = sqrt((sum((FS_P2$pol - FS_P2_moy)**2))/(length(FS_P2[,1])-1))
  Rep_P2 = 2*sqrt(2*Sr_P2)

  FS_P3_moy=mean(FS_P3$pol)
  Sr_P3 = sqrt((sum((FS_P3$pol - FS_P3_moy)**2))/(length(FS_P3[,1])-1))
  Rep_P3 = 2*sqrt(2*Sr_P3)

  Rep_all <- c(Rep_P1,Rep_P2,Rep_P3)
  Rep_all[Rep_all==0]=NA
  Rep_all=na.omit(Rep_all)

  Rep = mean(Rep_all)

  print(paste0("Repeatability for the ",i," sensor =",Rep))

 Rep_mean_all_sensors <- rbind(Rep_mean_all_sensors,Rep)   

}

Rep_sta = mean(Rep_mean_all_sensors)
print(paste0("Repeatability for ",pol," at the station ",sta_name[pp]," = ",round(Rep_sta)))


