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
#                    COMPARE SENSOR DATA AND MODEL                              #
#                Created 27/08/2019 updtaed 15/05/2020                          #
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
estim_period <- "112018"; estim_date <- "2018-11"; estim_date2 <- "201811" # estimation date
ndays <- 30 # number of days to be considered
unit="h" # unit of the comparison (time resolution h=hour)
estim_YYYY <- substr(estim_date,1,4); estim_MM <- substr(estim_date,6,7); estim_DD <- substr(estim_date,9,10)
estim_HH_start_list <- c("06","07","08","09","10","11","12","13","14","15","16","17","18")
estim_HH_end_list <- c("07","08","09","10","11","12","13","14","15","16","17","18","19")
all_stat_pol_model <- matrix(,nrow = length(estim_HH_start_list), ncol = 3)
all_stat_pol_data1 <- matrix(,nrow = length(estim_HH_start_list), ncol = 3)
all_stat_pol_data2 <- matrix(,nrow = length(estim_HH_start_list), ncol = 3)
R <- 6371e3 # medium earth radius in meters
X1=345000; X2=365000; Y1=6680000; Y2=6701000 # estimation domain limits in Lambert 93
all_obs_modelf <- c()
dist_max <- 5e9 # distance between points for comparison
CRS_L93="+init=epsg:2154" # Lambert 93 proj
CRS_WGS84="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # WGS84 proj

#####################################
#       READ RAW SENSOR DATA        #
#####################################

print("READ SENSOR DATA")

staFile <- paste0(indir,Msens_file)
stagrid <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
colnames(stagrid)[which(names(stagrid) == tolower(pol))] <- "pol"
stagrid <- subset(stagrid, select=-c(temperature, humidity,pressure))
nsensor <- unique(stagrid$id_sensor) # Sensors ID
stagrid$datetime = as.POSIXct(as.character(stagrid$datetime),tz='GMT')
tmp=stagrid$datetime
tmp2=as.character(tmp)
date=substr(tmp2,1,10)
hour=substr(tmp2,12,19)
hour_short=substr(tmp2,12,13)
hour_short_num=as.numeric(hour_short)
stagrid$hour_num <- hour_short_num
stagrid_ms_avg_hour <- aggregate(.~lat+lon+hour_num,data=stagrid,mean,na.action=na.omit) # hourly average

#####################################
#    READ CORRECTED SENSOR DATA     #
#####################################

print("READ CORRECTED SENSOR DATA")

# Init storage vector
cor_data_ms_all <- c()

for (dday in 1:ndays){ # loop over days to average for the entire period
    if (dday < 10){
        estim_DD <- paste0("0",toString(dday))}else{estim_DD <- toString(dday)}
    if (dday==4 | dday==11 | dday==18){
        cor_data_ms = cor_data_ms
    }else{
        file=paste("data_preproc_MS_bgdata_0.5_",estim_YYYY,estim_MM,estim_DD,".Rda",sep='')
        load(paste0(indir2,file)) # load data file
        cor_data_ms <- data_preproc
        cor_data_ms_all <- rbind(cor_data_ms,cor_data_ms_all)
    }
}

tmp=cor_data_ms_all$datetime
tmp2=as.character(tmp)
date=substr(tmp2,1,10)
hour=substr(tmp2,12,19)
hour_short=substr(tmp2,12,13)
hour_short_num=as.numeric(hour_short)
cor_data_ms_all$hour_num <- hour_short_num
cor_data_avg_hour <- aggregate(.~lat+lon+hour_num,data=cor_data_ms_all,mean,na.action=na.omit) # hourly average 
cor_data_avg_hour = subset(cor_data_avg_hour, select=-c(run))

#####################################
#       DO HOURLY COMPARISON        #
#####################################

print("START HOURLY COMPARISON BETWEEN DATA AND MODEL")

# Init variables
nbr_points_all <- c()

for (time in 1:length(estim_HH_start_list)){ # loop over time (hours of a day)

    estim_HH_start = estim_HH_start_list[time]
    estim_HH_end = estim_HH_end_list[time]
    estim_period = paste0(estim_MM,estim_YYYY,"_",estim_HH_start,"h",estim_HH_end,"h")
    print(paste0("TIME is ",estim_date," ", estim_HH_start,":00:00"))
    num_HH_start = as.numeric(estim_HH_start)
    num_HH_end = as.numeric(estim_HH_end)
    if (num_HH_end <= 9){ hhour=paste0("0",num_HH_end)}else{hhour=toString(num_HH_end)}

### Load model outputs ###

file <- paste0("ADMS_112018_",num_HH_end,"h_PM10_Nantes_7m.csv") # !!! MUST BE CHANGED BY THE USER !!!
data <- read.table(paste(indir,file,sep=''),header=T,sep=',',stringsAsFactors=F,dec='.')

spmodel=data # Convert to spatial data frame
names(spmodel)[names(spmodel) == "X"] <- "lon"
names(spmodel)[names(spmodel) == "Y"] <- "lat"
spmodel$Long <- spmodel$lon; spmodel$Lat <- spmodel$lat
coordinates(spmodel)=~Long+Lat
proj4string(spmodel)=CRS(CRS_L93) # define L93 as projection
spmodel <- spTransform(spmodel,CRS(CRS_WGS84)) # data transform to WGS84
spmodel <- subset(spmodel, select=-c(lat,lon))
tmp = spmodel@coords; dlon=tmp[,1]; dlat=tmp[,2]
spmodel$lon <- dlon
spmodel$lat <- dlat
class(spmodel); summary(spmodel)
names(spmodel)[names(spmodel) == pol] <- "pol"
spmodel=subset(spmodel,select=-c(num,NO2)) # !!! MUST BE UPDATED BY THE USER !!!
data_model = as.data.frame(spmodel)
names(data_model)[names(data_model) == "Long"] <- "lon"
names(data_model)[names(data_model) == "Lat"] <- "lat"
data_model = data_model[,c("lat","lon","pol")]

### Select corrected data ###

cor_data_sub_ms = subset(cor_data_avg_hour, hour_num >= num_HH_start & hour_num < num_HH_end) # subset data depending on the estimation period
cor_data_sub_ms_agg <- aggregate(.~lat+lon,data=cor_data_sub_ms,mean,na.action=na.pass) # average over lon/lat to avoid duplicated data

spdata_ms2 <- cor_data_sub_ms_agg # transform to spatial data frame 
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
spdata_ms2 = subset(spdata_ms2,select=-c(id_sensor,datetime,pm25,pol,hour_num)) # !!! MUST BE UPDATED BY THE USER !!!
names(spdata_ms2)[names(spdata_ms2) == "pol_bgcorr"] <- "pol"
data_ms2 = as.data.frame(spdata_ms2)
data_ms2 = subset(data_ms2,select=-c(Long,Lat))
data_ms2 = data_ms2[,c("lat","lon","pol")]

### Select raw sensor data ###

stagrid_ms_agg = subset(stagrid_ms_avg_hour, hour_num >= num_HH_start & hour_num < num_HH_end) # subset data depending on the estimation period
stagrid_sub_ms_agg <- aggregate(.~lat+lon,data=stagrid_ms_agg,mean,na.action=na.pass) # average over lon/lat to avoid duplicated data

stagrid_sub_ms_new <- stagrid_sub_ms_agg[1,] # select raw data depending on corrected data remaining after correction

for (mm in 1:length(cor_data_sub_ms_agg[,1])){
    idx=which(stagrid_sub_ms_agg$lat == cor_data_sub_ms_agg$lat[mm] & stagrid_sub_ms_agg$lon == cor_data_sub_ms_agg$lon[mm] )
    tmp=stagrid_sub_ms_agg[idx,]
    stagrid_sub_ms_new[mm,]=tmp
}

spdata_ms <- aggregate(.~lat+lon,data=stagrid_sub_ms_new,mean,na.action=na.pass) # mean over position to avoid duplicate data on a similar location ==> error for kriging
spdata_ms$Long <- spdata_ms$lon; spdata_ms$Lat <- spdata_ms$lat
coordinates(spdata_ms)=~Long+Lat
proj4string(spdata_ms)=CRS(CRS_WGS84)
spdata_ms=spTransform(spdata_ms,CRS(CRS_L93)) # data transform to L93
spdata_ms <- subset(spdata_ms, select=-c(lat,lon))
tmp = spdata_ms@coords; dlon=tmp[,1]; dlat=tmp[,2]
spdata_ms$lon <- dlon
spdata_ms$lat <- dlat
spdata_ms <- spdata_ms[spdata_ms$lon >= X1 & spdata_ms$lon <= X2 & spdata_ms$lat >= Y1 & spdata_ms$lat <= Y2,] # subset data to the grid limits
spdata_ms=spTransform(spdata_ms,CRS(CRS_WGS84)) # data transform to WGS84
tmp = spdata_ms@coords; dlon=tmp[,1]; dlat=tmp[,2]
spdata_ms$lon <- dlon
spdata_ms$lat <- dlat
class(spdata_ms); summary(spdata_ms)
spdata_ms = subset(spdata_ms,select=-c(id_sensor,datetime,pm25,hour_num)) # !!! MUST BE UPDATED BY THE USER !!!
data_ms = as.data.frame(spdata_ms)
data_ms = subset(data_ms,select=-c(Long,Lat))
data_ms = data_ms[,c("lat","lon","pol")]

# look at the number of data for representativity issue
nbr_points = length(data_ms[,1])
nbr_points_all = rbind(nbr_points_all,nbr_points)
 
### look for the closest model grid point to the data to be compared

model_tmp <- rep(0, length(data_ms[,1]))
for (ll in 1:length(spdata_ms2)){ # look for the closest model grid point to the data to be compared
    spdata_tmp = spdata_ms2[ll,]
    dist_vector = spDistsN1(spmodel,spdata_tmp) # calculate distance between the data point and the model grid points
    dist_vector = dist_vector * 1e3 # km to m
    dist_min=min(dist_vector)
    dist_min=min(dist_vector)
    idx=which.min(dist_vector) # look at the shortest distance
    tmp=spmodel[idx,]
    tmp2=as.data.frame(tmp)
    pol_tmp2=tmp2$pol
    model_tmp[ll]=pol_tmp2
}

obs_model <- cbind(data_ms$lon,data_ms$lat,data_ms$pol,model_tmp)
model_tmp <-subset(model_tmp, (!is.na(obs_model[,4])))
data_ms <-subset(data_ms, (!is.na(obs_model[,4])))
data_ms2 <-subset(data_ms2, (!is.na(obs_model[,4])))
obs_model <-subset(obs_model, (!is.na(obs_model[,4])))

### prepare data to plot with ggplot2 ###
LAT <- as.vector(t(rbind(data_ms$lat,data_ms$lat,data_ms$lat)))
LON <- as.vector(t(rbind(data_ms$lon,data_ms$lon,data_ms$lon)))
obs_model2 <- as.vector(t(rbind(data_ms$pol,data_ms2$pol,model_tmp)))
VAR <- as.vector(t(rbind(rep("Raw data",length(data_ms[,1])),rep("Corrected data",length(data_ms[,1])),rep("Model",length(data_ms[,1])))))
num_pts <- seq(from=1,to=length(data_ms[,1]))
PTS <-as.vector(t(rbind(num_pts,num_pts,num_pts)))
obs_modelf <- as.data.frame(cbind(PTS,LAT,LON,obs_model2))
obs_modelf$VAR <- as.factor(VAR)

all_obs_modelf <- rbind(all_obs_modelf,obs_modelf)

# Calculate average of all points and 95% confidence interval
n=length(data_ms[,1])

moy_pol_model = mean(model_tmp); moy_pol_data1 = mean(data_ms$pol); moy_pol_data2 = mean(data_ms2$pol)
sd_pol_model = sd(model_tmp); sd_pol_data1 = sd(data_ms$pol); sd_pol_data2 = sd(data_ms2$pol)

error_model <- qt(0.975,df=n-1)*sd_pol_model/sqrt(n)
error_data1 <- qt(0.975,df=n-1)*sd_pol_model/sqrt(n)
error_data2 <- qt(0.975,df=n-1)*sd_pol_model/sqrt(n)

upper_pol_model = moy_pol_model + error_model
upper_pol_data1 = moy_pol_data1 + error_data1
upper_pol_data2 = moy_pol_data2 + error_data2

lower_pol_model = moy_pol_model - error_model
lower_pol_data1 = moy_pol_data1 - error_data1
lower_pol_data2 = moy_pol_data2 - error_data2

all_stat_pol_model[time,1] = moy_pol_model
all_stat_pol_data1[time,1] = moy_pol_data1
all_stat_pol_data2[time,1] = moy_pol_data2

all_stat_pol_model[time,2] = upper_pol_model
all_stat_pol_data1[time,2] = upper_pol_data1
all_stat_pol_data2[time,2] = upper_pol_data2

all_stat_pol_model[time,3] = lower_pol_model
all_stat_pol_data1[time,3] = lower_pol_data1
all_stat_pol_data2[time,3] = lower_pol_data2

}

all_stat_pol <- as.data.frame(rbind(all_stat_pol_model,all_stat_pol_data1,all_stat_pol_data2))
names(all_stat_pol)[names(all_stat_pol) == "V1"] <- "Mean"
names(all_stat_pol)[names(all_stat_pol) == "V2"] <- "Upper"
names(all_stat_pol)[names(all_stat_pol) == "V3"] <- "Lower"
VAR1 <- rep("Model",length(all_stat_pol_model[,1]))
VAR2 <- rep("Raw data",length(all_stat_pol_model[,1]))
VAR3 <- rep("Corrected data",length(all_stat_pol_model[,1]))
VAR <- as.vector(t(rbind(VAR1,VAR2,VAR3)))
all_stat_pol$VAR <- as.factor(VAR)
ntime <- as.numeric(estim_HH_end_list)
TIME = as.vector(t(rbind(ntime,ntime,ntime)))
all_stat_pol$TIME <- TIME

# Do plot
source('multiplot.r')
cbPalette <- c("springgreen4","black","firebrick2")
cbPalette2 <- c("springgreen2","gray52","firebrick1")
cbPalette3 <- c("springgreen4","gray47","firebrick4")

png(paste0(dirout,"Figure2.png"),width=1000, height=500)
p1 <- ggplot(all_stat_pol,aes(x=TIME, y=Mean))+
geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=VAR), alpha=0.2) + 
geom_line(aes(color=VAR))+
scale_x_continuous(breaks = pretty(all_stat_pol$TIME, n = 12)) +
scale_fill_manual(values=cbPalette2)+
scale_colour_manual(values=cbPalette3)+
ggtitle("a)")+
xlab("Hour") + ylab(bquote(.(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p1
p2 <- ggplot(all_obs_modelf, aes(x=obs_model2, fill=VAR, color=VAR))+
geom_density(alpha=0.2)+
xlim(c(0, 100))+
scale_color_manual(values=alpha(cbPalette,0.4))+
scale_fill_manual(values=alpha(cbPalette,0.4))+
ggtitle("b)") +
xlab(bquote(.(pol) ~ (mu*g/m^3)))+
theme_bw()+
theme_minimal()+
theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_blank(),
        legend.spacing.x = unit(0.3, 'cm'),
        legend.position= "top")
p2
pf <- multiplot(p1,p2,cols=2,fontsize = 24)
dev.off()

