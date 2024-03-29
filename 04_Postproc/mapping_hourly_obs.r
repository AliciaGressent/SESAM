# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/04_Postproc/") # !!! MUST BE UPDATED BY THE USER !!!

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
#                     PLOT SENSOR DATA AT HOURLY RESOLUTION                     #
#                     Created 12/11/2018 updated 20/05/2020                     #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

# Define projections
CRS_L93=CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
CRS_WGS84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Situation
city="Nantes"
pol="PM10"

# Estimation time
estim_period="112018"
estim_date = "2018-11-29"
estim_date2 <- "29112018"
estim_YYYY <- substr(estim_date,1,4)
estim_MM <- substr(estim_date,6,7)
estim_DD <- substr(estim_date,9,10)
estim_HH_start_list <- c("06","07","08","09","10","11","12","13","14","15","16","17","18")
estim_HH_end_list <- c("07","08","09","10","11","12","13","14","15","16","17","18","19")
estim_HH_start_EN_list <- c("6","7","8","9","10","11","0","1","2","3","4","5","6")
estim_HH_end_EN_list <- c("7","8","9","10","11","0","1","2","3","4","5","6","7")
unit_list <- c("am","am","am","am","am","am","pm","pm","pm","pm","pm","pm","pm")

# Directory paths
indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/INPUTS/"          # Directory for input files
indir2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/"        # Directory2 for input files
file_drift <- "Drift_PM10_Nantes_7m.csv"                                                              # Drift file
file_ms <- paste("data_preproc_MS_bgdata_0.5_",estim_YYYY,estim_MM,estim_DD,".Rda",sep='')          # Mobile sensor data file
file_fs <- paste("data_preproc_FS_bgdata_0.5_",estim_YYYY,estim_MM,estim_DD,".Rda",sep='')          # Fixed sensor data file
outdir <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/04_Postproc/figs/" # Directory for output figures
outdir2 <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/"        # Directory for output files

# Parameters for grid definition
X1=345000; X2=365000; Y1=6680000; Y2=6701000    # Estimation domain limits
xmin=326361.53;                                 # Longitude min in Lambert 93
ymin=6675394.0;                                 # Latitude min in Lambert 93
res=7                                           # Grid Resolution in meters
ncols= 6619; nrows= 4466                        # Grid dimensions

# Shapefiles
roads_file <- "LOC_TRONCONS_ROUTIERS_NM.shp"
roads_shp  <- readShapeSpatial(paste0(indir,roads_file))
proj4string(roads_shp)=CRS_L93
e <- extent(344999,364998,6679996,6701004)
rc <- crop(roads_shp, e)

# Init for plot
png(paste0(outdir,"Figure8a.png"),width=1200, height=1200)
par(mar=c(5.1, 6.1, 4.1, 2.1))
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15), 4, 4, byrow = TRUE))

#####################################
#      CALCULATE VME AND PLOT       #          
#####################################

# Loop over hours
for (i in 1:length(estim_HH_start_list)){

        estim_HH_start = estim_HH_start_list[i]
        estim_HH_end = estim_HH_end_list[i]
        estim_date3 = paste0(estim_DD,estim_MM,estim_YYYY,"_",estim_HH_start,"h",estim_HH_end,"h")
        estim_HH_start_EN = estim_HH_start_EN_list[i]
        estim_HH_end_EN = estim_HH_end_EN_list[i]
        unit <- unit_list[i]
        file_EDK = paste0("EDK_grid_",city,"_",pol,"_",estim_date3,"_MS_FS.csv")

        # Select data for the period of estimation
        start_date=as.POSIXct(paste0(estim_date," ", estim_HH_start,":00:00"),origin = "1970-01-01", tz="GMT") # define period of estimation (1hour)
        end_date=as.POSIXct(paste0(estim_date," ", estim_HH_end,":00:00"),origin = "1970-01-01", tz="GMT")

        # Data from mobile sensors
        print("READ MOBILE SENSOR DATA")

        load(paste(indir2,file_ms,sep="")) # sensor data from preprocessing
        data_tmp1 <- data_preproc
        names(data_tmp1)[names(data_tmp1) == "pol"] <- "pol_ori"
        data_tmp1 = subset(data_tmp1, selec=-c(run,pm25)) # !!! MUST BE UPDATED BY THE USER !!! 
        data_tmp1 = na.omit(data_tmp1)
        data_tmp1 = subset(data_tmp1, datetime >= start_date & datetime <= end_date)
        if (length(data_tmp1[,1]) > 0){ # if data ms exist
            data_ms <- aggregate(.~lat+lon+id_sensor,data=data_tmp1,mean,na.action=na.pass)
            names(data_ms)[names(data_ms) == "pol_bgcorr"] <- "pol"
            names(data_ms)[names(data_ms) == "datetime"] <- "date"
        } # if data ms exist

        # Data from fixed sensors
        print("READ FIXED SENSOR DATA")

        load(paste(indir2,file_fs,sep="")) # sensor data from preprocessing
        data_tmp2 <- data_preproc
        names(data_tmp2)[names(data_tmp2) == "pol"] <- "pol_ori"
        data_tmp2 = subset(data_tmp2, selec=-c(run,pm25)) # !!! MUST BE UPDATED BY THE USER !!!
        data_tmp2=na.omit(data_tmp2)
        data_tmp2 = subset(data_tmp2, datetime >= start_date & datetime <= end_date)
        if (length(data_tmp2[,1]) > 0 ){ # if data fs exist
            data_fs <- aggregate(.~lat+lon+id_sensor,data=data_tmp2,mean,na.action=na.pass)
            names(data_fs)[names(data_fs) == "pol_bgcorr"] <- "pol"
            names(data_fs)[names(data_fs) == "datetime"] <- "date"
        } # if data fs exist

        # Concatenate data
        if (length(data_tmp1[,1]) > 0 & length(data_tmp2[,1]) > 0){
            data_tmp <- rbind(data_ms,data_fs); kriging_data <- 'MS_FS' # mobile and fixed sensor data
        }else if (length(data_tmp1[,1]) > 0 & length(data_tmp2[,1]) == 0){
            data_tmp <- rbind(data_ms); kriging_data <- 'MS' # only mobile sensor data
        }else if (length(data_tmp1[,1]) == 0 & length(data_tmp2[,1]) > 0){
            data_tmp <- rbind(data_fs); kriging_data <- 'FS' # only fixed sensor data
        }
       
        # Aggregate fixed and mobile sensor data and convert to spatial dataframe
        data <- aggregate(.~lat+lon,data=data_tmp,mean,na.action=na.pass) # mean over measurement position
        data$Long <- data$lon; data$Lat <- data$lat
        coordinates(data)=~Long+Lat
        proj4string(data)=CRS_WGS84
        data <- spTransform(data,CRS_L93) # data transform to Lambert 93
        data <- subset(data, select=-c(lat,lon))
        tmp = data@coords; dlon=tmp[,1]; dlat=tmp[,2]
        data$lon <- dlon
        data$lat <- dlat
        class(data); summary(data)
        data <- data[data$lon >= X1 & data$lon <= X2 & data$lat >= Y1 & data$lat <= Y2,] # subset data to the grid limits

        # Load estimate grid
        grid_EDK <-read.csv(paste0(indir2,file_EDK),header=TRUE,sep=",",skip=0)
        ras_pred_EDK=raster(list(x=sort(unique(coordinates(grid_EDK)[,1])),y=sort(unique(coordinates(grid_EDK)[,2])),z=matrix(grid_EDK$Pred_EDK,nrow=length(sort(unique(coordinates(grid_EDK)[,1]))),byrow=F)))

        ### PLOT

        # Color and scale for raster
        d=rbind(c(38,48,132),c(61,99,174),c(114,201,195),c(220,225,30),c(240,78,34),c(133,22,24))
        palette=colorRampPalette(rgb(d[,1],d[,2],d[,3],maxColorValue=255))(128)
        mycol=palette

        if (pol=="PM10"){ zl <- c(15,25)}
        ncol=100 # color and scale
        brks    <- seq(zl[1], zl[2], length.out = ncol+1)
        brkslab <- format(brks, scientific=FALSE, digits=2)
        indbrks <-  seq(1,length(brks), by = 15)
        mycol_ncol <- colorRampPalette(mycol)(ncol)

        data_plt <- data

        #Color and scale for data
        data_plt@data$pol[which(data_plt@data$pol <= zl[1])] <- zl[1]+0.1
        data_plt@data$pol[which(data_plt@data$pol >= zl[2])] <- zl[2]-0.1
        cut_pol <- cut(data_plt@data$pol,brks)
        col_pol <- mycol_ncol[cut_pol]

        # Plot
        plot(ras_pred_EDK, col="white", zlim=zl, breaks=brks, interpolate=FALSE, maxpixels=500000000,
            xaxt="n",yaxt="n", xaxs="i", yaxs="i",
            main=paste0(estim_HH_start_EN,'-',estim_HH_end_EN,unit), cex.main=2,
            legend=FALSE)
        plot(rc,add=TRUE,col="grey50",alpha=0.3)
        plot(data_plt,pch=21,col=col_pol,bg=col_pol,cex=2,add=TRUE)
}

dev.off()

png(paste0(outdir,"Figure8b.png"),width=1200, height=1200)
plot.new()
image.plot(legend.only=TRUE, horizontal = TRUE, zlim=zl,breaks=brks,
       col=mycol_ncol, legend.shrink=1,legend.width=2.5,legend.mar=50,
       legend.args=list(text=bquote(.(pol) ~ (mu*g/m^3)),cex=1.5,line=1,font=1,col='black'),
       axis.args=list(at=brks[indbrks],labels=format(brkslab[indbrks],digits=2),cex.axis=1.5,col.axis="black",col="black"))
dev.off()


