# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/04_Postproc/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(maptools);library(RColorBrewer)
library(fields);library(raster)
library(rgdal);library(gstat)
library(XML);library(akima)
library(rgeos);library(sp)
library(spacetime);library(chron)
library(plotrix); library(Rcpp)
library(RcppArmadillo); library(MASS)
library(ggplot2)

#################################################################################
#           COMPARE DATA FUSION TO OBSERVATIONS AND MODEL RESULTS               #
#                       Plot at each reference station                          #               
#                     Created 12/11/2018 updated 19/05/2020                     #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

print("INITIALIZATION")

# Define projections
CRS_L93=CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
CRS_WGS84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Situation
city="Nantes"
pol="PM10"
estim_period="112018"
estim_date1 = "2018-11-29" # estimation date
estim_date2 <- "29112018"
estim_YYYY <- substr(estim_date1,1,4) # extract year
estim_MM <- substr(estim_date1,6,7)   # extract month
estim_DD <- substr(estim_date1,9,10)  # extract day
estim_HH_start_list <- c("06","07","09","10","11","12","13","14","15","16","17","18") # starting hour
estim_HH_end_list <- c("07","08","10","11","12","13","14","15","16","17","18","19") # ending hour
time_start = as.POSIXct("2018-11-29 06:00:00",origin = "1970-01-01", tz="GMT") # start data selection
time_end = as.POSIXct("2018-11-29 19:00:00",origin = "1970-01-01", tz="GMT") # end data selection

# Directory paths and file names
indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/INPUTS/" # directory of inputs
indir2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/OUTPUTS/" # directory2 of inputs (mapping outputs) 
outdir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/OUTPUTS/" # directory of outputs
outdir2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/04_Postproc/figs/" # directory of plots
model_file <- "NA_ADMS-urban_nov2018v2.pst" # hourly output ADMS for the entire estimation period (November 2018)
ref_file <- "mesures_ref_qthourly_PM10.txt" # reference data file

# Other inputs
code_sta_pol <-c("188","238","107","239","140")
station_name <- c("Bouteillerie","Trentemoult","La Chauvinière","Les Couêts","Victor Hugo")
lat1 = 47.223253; lon1 = -1.539075 # Station N°1
lat2 = 47.194049; lon2 = -1.581454 # Station N°2
lat3 = 47.252721; lon3 = -1.574653 # Station N°3
lat4 = 47.185581; lon4 = -1.591332 # Station N°4
lat5 = 47.204264; lon5 = -1.552927 # Station N°5

#####################################
#        READ MODEL OUTPUTS         #          
#####################################

print("READ MODEL OUTPUTS")

# Read hourly model outputs
df_model <-read.table(paste0(indir,model_file),header=TRUE, sep=",",skip=0)
df_model_day <- df_model[which(df_model$Day ==333 & df_model$Hour > 5 & df_model$Hour < 20),] # select estimation date !!! MUST ME UPDATED BY THE USER !!!
df_pol_model = subset(df_model_day,select=-c(Day,Point,NO2,O3))
colnames(df_pol_model)[which(names(df_pol_model) == pol)] <- "pol"

#####################################
#          READ REF DATA            #          
#####################################

print("REFERENCE DATA")

# Read reference station measurements of the pollutant
refpol <- read.table(file=paste0(indir,ref_file), header=TRUE, sep="\t",skip=0)
refpol <- subset(refpol, select=-c(nom_station,polluant))
colnames(refpol)[which(names(refpol) == "niveau")] <- "pol"
refpol$date <- gsub("T", " ", refpol$date)
refpol$date = as.POSIXct(as.character(refpol$date),tz='GMT')
refpol=refpol[!duplicated(refpol),] # remove duplicated data
refpol <- na.omit(refpol) # remove NAN
nstations <- unique(refpol$code_station) # number of stations
lat_refpol <- rep(0,length(refpol[,1])); lon_refpol <- rep(0,length(refpol[,1]))
lon_refpol[which(refpol$code_station==code_sta_pol[1])]=lon1; lat_refpol[which(refpol$code_station==code_sta_pol[1])]=lat1;
lon_refpol[which(refpol$code_station==code_sta_pol[2])]=lon2; lat_refpol[which(refpol$code_station==code_sta_pol[2])]=lat2;
lon_refpol[which(refpol$code_station==code_sta_pol[3])]=lon3; lat_refpol[which(refpol$code_station==code_sta_pol[3])]=lat3;
lon_refpol[which(refpol$code_station==code_sta_pol[4])]=lon4; lat_refpol[which(refpol$code_station==code_sta_pol[4])]=lat4;
lon_refpol[which(refpol$code_station==code_sta_pol[5])]=lon5; lat_refpol[which(refpol$code_station==code_sta_pol[5])]=lat5;
refpol <- cbind(refpol,lat_refpol,lon_refpol)
names(refpol)[names(refpol) == "lat_refpol"] <- "lat"
names(refpol)[names(refpol) == "lon_refpol"] <- "lon"
df_ref_sub = subset(refpol, date >= time_start & date <= time_end)

# Hourly average
dtparts = t(as.data.frame(strsplit(as.character(df_ref_sub$date),' ')))
times = chron(dates=dtparts[,1],times=dtparts[,2], format=c('y-m-d','h:m:s'))
options(digits=12)
index_t = as.numeric(times) # number of days since 01/01/1970
df_ref_tmp <- df_ref_sub
df_ref_tmp$index_t <- index_t
averageT1 = 60*6.94e-4 # minute convert in days (here 1 hour average)
timevec_out1 <- seq(from=min(unique(df_ref_tmp$index_t)/averageT1),to=max(unique(df_ref_tmp$index_t)/averageT1),by=1)*averageT1
df_ref_tmp$intervalT1 <- cut(df_ref_tmp$index_t,breaks=timevec_out1)
df_ref_hours <- aggregate(. ~ lon+lat+intervalT1,data=df_ref_tmp,mean,na.action=na.pass)
df_ref_hours$date = as.POSIXct(df_ref_hours$date,origin = "1970-01-01", tz="GMT")
df_ref_hours = subset(df_ref_hours,select=-c(intervalT1,index_t))
df_ref_hours$Long <- df_ref_hours$lon; df_ref_hours$Lat <- df_ref_hours$lat
coordinates(df_ref_hours)=~Long+Lat
proj4string(df_ref_hours)=CRS_WGS84
df_ref_hours <- spTransform(df_ref_hours,CRS_L93) # data transform to Lambert 93
df_ref_hours = as.data.frame(df_ref_hours) # convert to a dataframe

# Daily average
df_ref_day <- aggregate(.~lat+lon,data=df_ref_sub,mean,na.action=na.pass)

#####################################
#   SELECT DATA AND MODE HOURLY     #          
#####################################

print("START SELECT PRED, MODEL AND DATA AT HOURLY RESOLUTION")

# Init storage vectors
pol_pred_obs <- c()
pol_model_obs <- c()

# Loop over hours
for (i in 1:length(estim_HH_start_list)){

    # Define data fusion estimation file
    estim_HH_start = estim_HH_start_list[i]
    estim_HH_end = estim_HH_end_list[i]
    estim_date3 = paste0(estim_DD,estim_MM,estim_YYYY,"_",estim_HH_start,"h",estim_HH_end,"h")
    data_type="MS_FS"
    pred_file = paste0("EDK_grid_Nantes_",pol,"_",estim_date3,"_",data_type,".csv")
    print(pred_file)
    len_file <- nchar(pred_file)
    kriging_data <- substr(pred_file, 38, 42)
    estim_date <- substr(pred_file, 22, 36)

    # Select reference obs - hour
    hour_start = as.POSIXct(paste0(estim_date1," ",substr(pred_file, 31, 32),":00:00"),origin = "1970-01-01", tz="GMT")
    hour_end = as.POSIXct(paste0(estim_date1," ",substr(pred_file, 34, 35),":00:00"),origin = "1970-01-01", tz="GMT")
    df_ref_t = subset(df_ref_hours, date >= hour_start & date <= hour_end)
        
    # Subset model outputs for time period
    hour = substr(pred_file, 34, 35)
    if (hour == "06" | hour == "07" | hour == "08" | hour == "09"){ hour = substr(hour, 2, 2)}else{hour = hour}
    df_pol_model_h <- df_pol_model[which(df_pol_model$Hour == hour),] # select hour
        
    # Load estimation grids
    grid_EDK <-read.csv(paste0(indir2,pred_file),header=TRUE,sep=",",skip=0)
        
    # Mean grid and create a raster
    if (i == 1){
        grid_Pred_EDK <- grid_EDK
        grid_Pred_EDK = subset(grid_Pred_EDK,select=-c(Drift,Long,Lat,optional))
    }else{
        grid_Pred_EDK$Pred_EDK = grid_Pred_EDK$Pred_EDK + grid_EDK$Pred_EDK
        grid_Pred_EDK$StDev_EDK = grid_Pred_EDK$StDev_EDK + grid_EDK$StDev_EDK
    }

     # Select EDK estimation for the pollutant and model outputs at station locations
     pol_pred <- rep(0, length(df_ref_t[,1]))
     pol_model <- rep(0, length(df_ref_t[,1]))

     for (j in 1:length(df_ref_t[,1])){
         pol_pred[j]=grid_EDK[which.min(abs(grid_EDK$lon - df_ref_t[j,]$Lon) + abs(grid_EDK$lat - df_ref_t[j,]$Lat)),6]
         pol_model[j]=df_pol_model_h[which.min(abs(df_pol_model_h$x - df_ref_t[j,]$Lon) + abs(df_pol_model_h$y - df_ref_t[j,]$Lat)),4]
     }

     # Concatenate
     pol_pred_obs_tmp = cbind(df_ref_t$code_station,df_ref_t$lon,df_ref_t$lat,df_ref_t$date,df_ref_t$pol,pol_pred)
     pol_pred_obs = rbind(pol_pred_obs,pol_pred_obs_tmp)
     pol_model_obs_tmp = cbind(df_ref_t$code_station,df_ref_t$lon,df_ref_t$lat,df_ref_t$date,df_ref_t$pol,pol_model)
     pol_model_obs = rbind(pol_model_obs,pol_model_obs_tmp)

}

grid_Pred_EDK$Pred_EDK = grid_Pred_EDK$Pred_EDK/i # average pred over the estimation period
grid_Pred_EDK$StDev_EDK = grid_Pred_EDK$StDev_EDK/i # average error over the estimation period

write.table(grid_Pred_EDK,file=paste(outdir,"EDK_AVG_grid_",city,"_",pol,"_",estim_DD,estim_MM,estim_YYYY,".csv",sep=''),row.names=F,col.names=T,sep=',') # save estimation grid and values

ras_pred_EDK=raster(list(x=sort(unique(coordinates(grid_Pred_EDK)[,1])),y=sort(unique(coordinates(grid_Pred_EDK)[,2])),z=matrix(grid_Pred_EDK$Pred_EDK,nrow=length(sort(unique(coordinates(grid_Pred_EDK)[,1]))),byrow=F))) # create a raster

#####################################
#              DO PLOT              #          
#####################################

### Prepare data for plot

pol_pred_obs = as.data.frame(pol_pred_obs)
names(pol_pred_obs)[1] <-"code_station"; names(pol_pred_obs)[2] <-"lon"; names(pol_pred_obs)[3] <-"lat"
names(pol_pred_obs)[4] <-"date"; names(pol_pred_obs)[5] <-"pol_obs"; names(pol_pred_obs)[6] <-"pol_pred"
pol_pred_obs$date <- as.POSIXct(pol_pred_obs$date,origin = "1970-01-01", tz="GMT")
pol_pred_obs$diff <- pol_pred_obs$pol_obs - pol_pred_obs$pol_pred
pol_pred_obs$code_station[pol_pred_obs$code_station==code_sta_pol[1]]<-station_name[1]
pol_pred_obs$code_station[pol_pred_obs$code_station==code_sta_pol[2]]<-station_name[2]
pol_pred_obs$code_station[pol_pred_obs$code_station==code_sta_pol[3]]<-station_name[3]
pol_pred_obs$code_station[pol_pred_obs$code_station==code_sta_pol[4]]<-station_name[4]
pol_pred_obs$code_station[pol_pred_obs$code_station==code_sta_pol[5]]<-station_name[5]
pol_model_obs = as.data.frame(pol_model_obs)
names(pol_model_obs)[1] <-"code_station"; names(pol_model_obs)[2] <-"lon"; names(pol_model_obs)[3] <-"lat"; 
names(pol_model_obs)[4] <-"date"; names(pol_model_obs)[5] <-"pol_obs"; names(pol_model_obs)[6] <-"pol_model"
pol_model_obs$date <- as.POSIXct(pol_model_obs$date,origin = "1970-01-01", tz="GMT")
pol_model_obs$diff <- pol_model_obs$pol_obs - pol_model_obs$pol_model
pol_model_obs$code_station[pol_model_obs$code_station==code_sta_pol[1]]<-station_name[1]
pol_model_obs$code_station[pol_model_obs$code_station==code_sta_pol[2]]<-station_name[2]
pol_model_obs$code_station[pol_model_obs$code_station==code_sta_pol[3]]<-station_name[3]
pol_model_obs$code_station[pol_model_obs$code_station==code_sta_pol[4]]<-station_name[4]
pol_model_obs$code_station[pol_model_obs$code_station==code_sta_pol[5]]<-station_name[5]

### Do plot

# Init variables
var <- c("pol_sta1","pol_sta2","pol_sta3","pol_sta4","pol_sta5")
var_model <- c("pol_sta1_model","pol_sta2_model","pol_sta3_model","pol_sta4_model","pol_sta5_model")

### Hourly comparison at each station between data fusion, observations and model

for (n in 1:length(nstations)){
    ncode <- code_sta_pol[n]
    nname <- station_name[n]
    pol_sta <- pol_pred_obs[which(pol_pred_obs$code_station == nname),]
    pol_stab <- pol_model_obs[which(pol_model_obs$code_station == nname),]
    assign(var_model[n], pol_stab)
    assign(var[n], pol_sta)
}

png(filename=paste0(outdir2,"Figure1.png"), width=800, height=800, type="cairo",bg = "white")
par(mfrow=c(3,2), mar=c(5,5,5,3)+0.0)
lower_lim=0; upper_lim=40
plot(pol_sta1$date,pol_sta1$pol_obs, pch=19, cex=1.5,
        xlab="Date", cex.lab=2,
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim),panel.first = grid())
title(station_name[1],cex.main=2,adj=0,line=0.5)
points(pol_sta1$date,pol_sta1$pol_pred,col="Red",pch=19, cex=1.5)
points(pol_sta1_model$date,pol_sta1_model$pol_model,col="dodgerblue3",pch=19, cex=1.5)
legend("topleft", horiz=TRUE, legend=c("Ref. observations", "Fusion","Model"),
            col=c("black","red","dodgerblue3"), pch=c(19,19,19),cex=1.5,box.lty=0)

plot(pol_sta2$date,pol_sta2$pol_obs, pch=19, cex=1.5,
        xlab="Date", cex.lab=2,
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim),panel.first = grid())
title(station_name[2],cex.main=2,adj=0,line=0.5)
points(pol_sta2$date,pol_sta2$pol_pred,col="Red",pch=19, cex=1.5)
points(pol_sta2_model$date,pol_sta2_model$pol_model,col="dodgerblue3",pch=19, cex=1.5)
legend("topleft", horiz=TRUE, legend=c("Ref. observations", "Fusion","Model"),
            col=c("black","red","dodgerblue3"), pch=c(19,19,19),cex=1.5,box.lty=0)

plot(pol_sta3$date,pol_sta3$pol_obs, pch=19, cex=1.5,
        xlab="Date", cex.lab=2,
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim),panel.first = grid())
title(station_name[3],cex.main=2,adj=0,line=0.5)
points(pol_sta3$date,pol_sta3$pol_pred,col="Red",pch=19,cex=1.5)
points(pol_sta3_model$date,pol_sta3_model$pol_model,col="dodgerblue3",pch=19, cex=1.5)
legend("topleft", horiz=TRUE, legend=c("Ref. observations", "Fusion","Model"),
            col=c("black","red","dodgerblue3"), pch=c(19,19,19),cex=1.5,box.lty=0)

plot(pol_sta4$date,pol_sta4$pol_obs, pch=19,cex=1.5,
        xlab="Date", cex.lab=2,
        ylab=bquote(LCS ~ .(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim),panel.first = grid())
title(station_name[4],cex.main=2,adj=0,line=0.5)
points(pol_sta4$date,pol_sta4$pol_pred,col="Red",pch=19,cex=1.5)
points(pol_sta4_model$date,pol_sta4_model$pol_model,col="dodgerblue3",pch=19, cex=1.5)
legend("topleft", horiz=TRUE, legend=c("Ref. observations", "Fusion","Model"),
            col=c("black","red","dodgerblue3"), pch=c(19,19,19),cex=1.5,box.lty=0)

plot(pol_sta5$date,pol_sta5$pol_obs, pch=19, cex=1.5,
        xlab="Date", cex.lab=2,
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim),panel.first = grid())
title(station_name[5],cex.main=2,adj=0,line=0.5)
points(pol_sta5$date,pol_sta5$pol_pred,col="Red",pch=19, cex=1.5)
points(pol_sta5_model$date,pol_sta5_model$pol_model,col="dodgerblue3",pch=19, cex=1.5)
legend("topleft", horiz=TRUE, legend=c("Ref. observations", "Fusion","Model"),
            col=c("black","red","dodgerblue3"), pch=c(19,19,19),cex=1.5,box.lty=0)
dev.off()

### Daily comparison at each station between data fusion, observations and model

for (n in 1:length(nstations)){
    ncode <- code_sta_pol[n]
    nname <- station_name[n]

    # Data fusion estimation
    pol_sta_pred <- pol_pred_obs[which(pol_pred_obs$code_station == nname),] # select pred for the nieme station
    pol_sta_pred_mean <- mean(pol_sta_pred$pol_pred) # mean over the period
    pol_sta_pred_sd <- sd(pol_sta_pred$pol_pred) # stdev over the period

    # Model outputs
    pol_sta_model <- pol_model_obs[which(pol_model_obs$code_station == nname),] # select model for the nieme station
    pol_sta_model_mean <- mean(pol_sta_model$pol_model) # mean over the period
    pol_sta_model_sd <- sd(pol_sta_model$pol_model) # stdev over the period

    # Observations at monitoring stations
    pol_sta_obs_mean <- mean(pol_sta_model$pol_obs) # mean over the period
    pol_sta_obs_sd <- sd(pol_sta_model$pol_obs) # stdev over the period

    pol_sta_mean <- cbind(pol_sta_obs_mean,pol_sta_pred_mean,pol_sta_model_mean)
    pol_sta_sd <- cbind(pol_sta_obs_sd,pol_sta_pred_sd,pol_sta_model_sd)

    # Concatenate
    if (n==1){
        pol_sta_mean_all <- pol_sta_mean
        pol_sta_sd_all <- pol_sta_sd
    }else{
        pol_sta_mean_all = rbind(pol_sta_mean_all,pol_sta_mean)
        pol_sta_sd_all = rbind(pol_sta_sd_all,pol_sta_sd)
    }

}

pol_sta_sd_all = as.data.frame(pol_sta_sd_all)
pol_sta_sd_all = cbind(pol_sta_sd_all,station_name)
names(pol_sta_sd_all)[1] <-"Ref"
names(pol_sta_sd_all)[2] <-"Estimate"
names(pol_sta_sd_all)[3] <-"Model"
names(pol_sta_sd_all)[4] <-"Station"

pol_sta_mean_all = as.data.frame(pol_sta_mean_all)
pol_sta_mean_all = cbind(pol_sta_mean_all,station_name)
names(pol_sta_mean_all)[1] <-"Ref"
names(pol_sta_mean_all)[2] <-"Estimate"
names(pol_sta_mean_all)[3] <-"Model"
names(pol_sta_mean_all)[4] <-"Station"

png(filename=paste0(outdir2,"Figure2.png"), width=800, height=800, type="cairo",bg = "white")
par(mfrow=c(3,2), mar=c(5,5,5,3)+0.0)
lower_lim=10; upper_lim=40

mean_conc=as.numeric(pol_sta_mean_all[1,1:3]); stdev_conc=as.numeric(pol_sta_sd_all[1,1:3])
plot(mean_conc, pch=19, cex=1.5, col=c("black",'red',"dodgerblue3"), xaxt = 'n',
        xlab="",
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim), xlim=c(0.5,3.5),panel.first = grid())
arrows(c(1,2,3), mean_conc-stdev_conc, c(1,2,3), mean_conc+stdev_conc, length=0.05, angle=90, code=3,col=c("black",'red',"dodgerblue3"))
axis(1, cex.axis=2, at=1:3, labels=c("Ref. observations","Fusion","Model"))
title(station_name[1],cex.main=2,adj=0,line=0.5)

mean_conc=as.numeric(pol_sta_mean_all[2,1:3]); stdev_conc=as.numeric(pol_sta_sd_all[2,1:3])
plot(mean_conc, pch=19, cex=1.5, col=c("black",'red',"dodgerblue3"), xaxt = 'n',
        xlab="",
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim), xlim=c(0.5,3.5),panel.first = grid())
arrows(c(1,2,3), mean_conc-stdev_conc, c(1,2,3), mean_conc+stdev_conc, length=0.05, angle=90, code=3,col=c("black",'red',"dodgerblue3"))
axis(1, cex.axis=2, at=1:3, labels=c("Ref. observations","Fusion","Model"))
title(station_name[2],cex.main=2,adj=0,line=0.5)

mean_conc=as.numeric(pol_sta_mean_all[3,1:3]); stdev_conc=as.numeric(pol_sta_sd_all[3,1:3])
plot(mean_conc, pch=19, cex=1.5, col=c("black",'red',"dodgerblue3"), xaxt = 'n',
        xlab="",
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim), xlim=c(0.5,3.5),panel.first = grid())
arrows(c(1,2,3), mean_conc-stdev_conc, c(1,2,3), mean_conc+stdev_conc, length=0.05, angle=90, code=3,col=c("black",'red',"dodgerblue3"))
axis(1, cex.axis=2, at=1:3, labels=c("Ref. observations","Fusion","Model"))
title(station_name[3],cex.main=2,adj=0,line=0.5)

mean_conc=as.numeric(pol_sta_mean_all[4,1:3]); stdev_conc=as.numeric(pol_sta_sd_all[4,1:3])
plot(mean_conc, pch=19, cex=1.5, col=c("black",'red',"dodgerblue3"), xaxt = 'n',
        xlab="",
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim), xlim=c(0.5,3.5),panel.first = grid())
arrows(c(1,2,3), mean_conc-stdev_conc, c(1,2,3), mean_conc+stdev_conc, length=0.05, angle=90, code=3,col=c("black",'red',"dodgerblue3"))
axis(1, cex.axis=2, at=1:3, labels=c("Ref. observations","Fusion","Model"))
title(station_name[4],cex.main=2,adj=0,line=0.5)

mean_conc=as.numeric(pol_sta_mean_all[5,1:3]); stdev_conc=as.numeric(pol_sta_sd_all[5,1:3])
plot(mean_conc, pch=19, cex=1.5, col=c("black",'red',"dodgerblue3"), xaxt = 'n',
        xlab="",
        ylab=bquote(.(pol) ~ (mu*g/m^3)),cex.lab=2,cex.axis = 2,
        ylim=c(lower_lim,upper_lim), xlim=c(0.5,3.5),panel.first = grid())
arrows(c(1,2,3), mean_conc-stdev_conc, c(1,2,3), mean_conc+stdev_conc, length=0.05, angle=90, code=3,col=c("black",'red',"dodgerblue3"))
axis(1, cex.axis=2, at=1:3, labels=c("Ref. observations","Fusion","Model"))
title(station_name[5],cex.main=2,adj=0,line=0.5)

dev.off()

### Map of the daily average

# Color and scale for raster plot
d=rbind(c(38,48,132),c(61,99,174),c(114,201,195),c(220,225,30),c(240,78,34),c(133,22,24))
palette=colorRampPalette(rgb(d[,1],d[,2],d[,3],maxColorValue=255))(128)
mycol=palette

if (pol=="PM10"){ zl <- c(15,25)}
ncol=100 # number of colors
brks    <- seq(zl[1], zl[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)
ras_pred_EDK[which(ras_pred_EDK@data@values < zl[1]) ] <- zl[1]
ras_pred_EDK[which(ras_pred_EDK@data@values > zl[2]) ] <- zl[2]

data_ref= aggregate(.~lon+lat, data=pol_pred_obs, mean,na.action=na.omit)
data_ref$Long <- data_ref$lon; data_ref$Lat <- data_ref$lat
coordinates(data_ref)=~Long+Lat
proj4string(data_ref)=CRS_WGS84
data_ref <- spTransform(data_ref,CRS_L93) # data transform to Lambert 93
data_ref <- subset(data_ref, select=-c(lat,lon))
tmp = data_ref@coords; dlon=tmp[,1]; dlat=tmp[,2]
data_ref$lon <- dlon
data_ref$lat <- dlat
class(data_ref); summary(data_ref)

png(filename=paste0(outdir2,"Figure3.png"),width=1000, height=800)
par(mar=c(4,6,4,7)) # margin bot left top right (need space for the 
plot(ras_pred_EDK, col=mycol_ncol, zlim=zl, breaks=brks, interpolate=TRUE,
     xaxt="n",yaxt="n", xaxs="i", yaxs="i",
     main="c) Fused map", cex.main=3,
     #xlab="Longitude",ylab="Latitude", cex.lab=1.5, 
     cex.axis=2.5,
     legend.width = 1.5, legend.shrink=0.75,
     #legend.args=list(text="(ug/m3)",cex=2.5,line=1,font=1),
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=2.5,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=2.5))
#plot(rc,add=TRUE,col="grey50",alpha=0.7)
#plot(data_ref,pch=17,col='black',add=T,cex=1.5)
#plot(data_plt,pch=21,col=alpha('black', 0.4),bg=col_pm10,add=T,cex=2.5)
dev.off()



