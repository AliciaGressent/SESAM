# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/02_Preproc/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(raster)
library(maptools)
library(RColorBrewer)
library(fields)
library(rgdal)
library(ggplot2)
library(chron)
library(geosphere)

#################################################################################
#                CALCULATE CORRELATION BETWEEN DATA AND DRIFT                   #
#                     Created 17/10/2019 updated 15/05/2020                     #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

print("INITIALIZATION")

indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/INPUTS/" # path for input directory
indir2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/" # path for input directory2
dirout <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/02_Preproc/figs/" # path for output directory plot
dirout2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/" # path for output directory files 
pol <- "PM10" # pollutant
Msens_file <- "dataout_mobile_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
Fsens_file <- "dataout_fixe_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
drift_file <- "Drift_PM10_Nantes_7m.csv" # ADMS-urban outputs and interpolated: 2016 annual mean
estim_period <- "112018"; estim_date = "2018-11-29"; estim_date2="20181129"; # estimation date
estim_YYYY <- substr(estim_date,1,4); estim_MM <- substr(estim_date,6,7); estim_DD <- substr(estim_date,9,10)
estim_HH_start_list <- c("07","08","09","10","11","12","13","14","15","16","17","18")
estim_HH_end_list <- c("08","09","10","11","12","13","14","15","16","17","18","19")
all_stat_pol_model <- matrix(,nrow = length(estim_HH_start_list), ncol = 3)
all_stat_pol_data <- matrix(,nrow = length(estim_HH_start_list), ncol = 3)
R <- 6371e3 # Medium earth radius in meters
X1=345000; X2=365000; Y1=6680000; Y2=6701000 # Domain limits in L93
all_obs_modelf <- c()
dist_max=0
unit="h" # unit of the comparison (time resolution h=hour)
CRS_L93="+init=epsg:2154" # Lambert 93 proj
CRS_WGS84="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # WGS84 proj

#####################################
#    READ CORRECTED SENSOR DATA     #
#####################################

print("READ CORRECTED SENSOR DATA")

file=paste0("data_preproc_MS_bgdata_0.5_",estim_YYYY,estim_MM,estim_DD,".Rda")
load(paste0(indir2,file)) 
cor_data_ms <- data_preproc

file=paste0("data_preproc_FS_bgdata_0.5_",estim_YYYY,estim_MM,estim_DD,".Rda") 
load(paste0(indir2,file))
cor_data_fs <- data_preproc

#####################################
#            LOAD DRIFT             #
#####################################

print("READ DRIFT")

var_aux <-read.csv(paste0(indir2,drift_file),header=TRUE,sep=",",skip=0)

# Subset auxiliary variable dataset to the data "domain"
names(var_aux)[names(var_aux) == "X"] <- "lon"
names(var_aux)[names(var_aux) == "Y"] <- "lat"
var_auxb <- var_aux
var_auxc <- var_auxb
var_auxc = subset(var_auxb, lon >= X1 & lon <= X2 & lat >= Y1 & lat <= Y2) # Subset full grid to the data geographical limits
spmodel <- var_auxc
spmodel$Long <- spmodel$lon; spmodel$Lat <- spmodel$lat
coordinates(spmodel)=~Long+Lat
proj4string(spmodel)=CRS(CRS_L93)
spmodel <- spTransform(spmodel,CRS(CRS_WGS84)) # data transform to WGS84
spmodel <- subset(spmodel, select=-c(lat,lon))
tmp = spmodel@coords; dlon=tmp[,1]; dlat=tmp[,2]
spmodel$lon <- dlon
spmodel$lat <- dlat
class(spmodel); summary(spmodel)
names(spmodel)[names(spmodel) == "drift_pol"] <- "pol"
data_model = as.data.frame(spmodel)
names(data_model)[names(data_model) == "Long"] <- "lon"
names(data_model)[names(data_model) == "Lat"] <- "lat"
data_model = data_model[,c("lat","lon","pol")]

#####################################
#     LOOK AT THE CORRELATION       #
#####################################

# Init
corr_all <- c()

# Loop over time (houris of a day)
for (time in 1:length(estim_HH_start_list)){

    estim_HH_start = estim_HH_start_list[time]
    estim_HH_end = estim_HH_end_list[time]
    estim_period = paste0(estim_DD,estim_MM,estim_YYYY,"_",estim_HH_start,"h",estim_HH_end,"h")
    print(paste0("TIME is ",estim_date," ", estim_HH_start,":00:00"))
    num_HH_start = as.numeric(estim_HH_start)
    num_HH_end = as.numeric(estim_HH_end)
    if (num_HH_end <= 9){ hhour=paste0("0",num_HH_end)}else{hhour=toString(num_HH_end)}
    datestart=as.POSIXct(paste0(estim_date," ", estim_HH_start,":00:00"),origin = "1970-01-01", tz="GMT") # define period of estimation (1hour)
    datefin=as.POSIXct(paste0(estim_date," ", estim_HH_end,":00:00"),origin = "1970-01-01", tz="GMT")

### Select corrected sensor data ###

    ## Mobile sensor data ## 
    cor_data_sub_ms = subset(cor_data_ms, datetime >= datestart & datetime <= datefin) # subset data depending on the estimation period
    cor_data_sub_ms_agg <- aggregate(.~lat+lon,data=cor_data_sub_ms,mean,na.action=na.pass) # mean over position to avoid duplicate data on a similar location ==> error for kriging

    spdata_ms2 <- cor_data_sub_ms_agg
    spdata_ms2$Long <- spdata_ms2$lon; spdata_ms2$Lat <- spdata_ms2$lat
    coordinates(spdata_ms2)=~Long+Lat
    proj4string(spdata_ms2)=CRS(CRS_WGS84)
    spdata_ms2=spTransform(spdata_ms2,CRS(CRS_L93)) # data transform to L93
    spdata_ms2 <- subset(spdata_ms2, select=-c(lat,lon))
    tmp = spdata_ms2@coords; dlon=tmp[,1]; dlat=tmp[,2]
    spdata_ms2$lon <- dlon
    spdata_ms2$lat <- dlat
    spdata_ms2 <- spdata_ms2[spdata_ms2$lon >= X1 & spdata_ms2$lon <= X2 & spdata_ms2$lat >= Y1 & spdata_ms2$lat <= Y2,] # subset data to the grid limits
    spdata_ms2=spTransform(spdata_ms2,CRS(CRS_WGS84)) # data transform to WGS84
    tmp = spdata_ms2@coords; dlon=tmp[,1]; dlat=tmp[,2]
    spdata_ms2$lon <- dlon
    spdata_ms2$lat <- dlat
    class(spdata_ms2); summary(spdata_ms2)
    spdata_ms2 = subset(spdata_ms2,select=-c(id_sensor,datetime,pm25,pol,run)) # !!! MUST BE UPDATED BY THE USER !!!
    names(spdata_ms2)[names(spdata_ms2) == "pol_bgcorr"] <- "pol"
    data_ms2 = as.data.frame(spdata_ms2)
    data_ms2 = subset(data_ms2,select=-c(Long,Lat))
    data_ms2 = data_ms2[,c("lat","lon","pol")]

    ## Fixed sensor data ## 
    cor_data_sub_fs = subset(cor_data_fs, datetime >= datestart & datetime <= datefin) # subset data depending on the estimation period
    cor_data_sub_fs_agg <- aggregate(.~lat+lon,data=cor_data_sub_fs,mean,na.action=na.pass) # mean over position to avoid duplicate data on a similar location ==> error for kriging

    spdata_fs2 <- cor_data_sub_fs_agg
    spdata_fs2$Long <- spdata_fs2$lon; spdata_fs2$Lat <- spdata_fs2$lat
    coordinates(spdata_fs2)=~Long+Lat
    proj4string(spdata_fs2)=CRS(CRS_WGS84)
    spdata_fs2=spTransform(spdata_fs2,CRS(CRS_L93)) # data transform to L93
    spdata_fs2 <- subset(spdata_fs2, select=-c(lat,lon))
    tmp = spdata_fs2@coords; dlon=tmp[,1]; dlat=tmp[,2]
    spdata_fs2$lon <- dlon
    spdata_fs2$lat <- dlat
    spdata_fs2 <- spdata_fs2[spdata_fs2$lon >= X1 & spdata_fs2$lon <= X2 & spdata_fs2$lat >= Y1 & spdata_fs2$lat <= Y2,] # subset data to the grid limits
    spdata_fs2=spTransform(spdata_fs2,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # data transform to WGS84
    tmp = spdata_fs2@coords; dlon=tmp[,1]; dlat=tmp[,2]
    spdata_fs2$lon <- dlon
    spdata_fs2$lat <- dlat
    class(spdata_fs2); summary(spdata_fs2)
    spdata_fs2 = subset(spdata_fs2,select=-c(id_sensor,datetime,pm25,pol)) # !!! MUST BE UPDATED BY THE USER !!!
    names(spdata_fs2)[names(spdata_fs2) == "pol_bgcorr"] <- "pol"
    data_fs2 = as.data.frame(spdata_fs2)
    data_fs2 = subset(data_fs2,select=-c(Long,Lat))
    data_fs2 = data_fs2[,c("lat","lon","pol")]

    data_all <- rbind(data_ms2,data_fs2) # concatenate mobile and fixed sensor data
    data <- aggregate(.~lat+lon,data=data_all,mean,na.action=na.pass)

    spdata <- data # transform to spatial dataframe
    spdata$Long <- spdata$lon; spdata$Lat <- spdata$lat
    coordinates(spdata)=~Long+Lat
    proj4string(spdata)=CRS("+init=epsg:4326")
    spdata <- subset(spdata, select=-c(lat,lon))
    tmp = spdata@coords; dlon=tmp[,1]; dlat=tmp[,2]
    spdata$lon <- dlon
    spdata$lat <- dlat

### Select model points and caluclate correlation ###

    pol_model <- rep(0, length(data[,1]))

    for (ll in 1:length(spdata)){
        spdata_tmp = spdata[ll,]
        dist_vector = spDistsN1(spmodel,spdata_tmp,longlat = TRUE)
        dist_vector = dist_vector * 1e3 # km to m
        dist_min=min(dist_vector)        
        idx=which(dist_vector == dist_min)
        tmp=spmodel[idx,]
        tmp2=as.data.frame(tmp)
        pol_tmp2=mean(tmp2$pol) 
        pol_model[ll]=pol_tmp2
    }

    corr=round(cor(spdata$pol,pol_model),2)
    print(paste0("Correlation is ", corr))
    corr_all = rbind(corr_all,corr)

} # Loop over time (hours of a day)

save(corr_all,file=paste0(dirout,"corr_",dist_max,"m_all_",estim_date2,".Rda"))
