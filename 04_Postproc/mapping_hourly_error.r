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
#                  PLOT DATA FUSION ERROR AT HOURLY RESOLUTION                  #
#                     Created 12/11/2018 updated 25/05/2020                     #
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
outdir <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/04_Postproc/figs/" # Directory for output figures
outdir2 <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/"        # Directory for output files

# Parameters for grid definition
X1=345000; X2=365000; Y1=6680000; Y2=6701000    # Estimation domain limits
xmin=326361.53;                                 # Longitude min in Lambert 93
ymin=6675394.0;                                 # Latitude min in Lambert 93
res=7                                           # Grid Resolution in meters
ncols= 6619; nrows= 4466                        # Grid dimensions

# Shapefile
roads_file <- "LOC_TRONCONS_ROUTIERS_NM.shp"
roads_shp  <- readShapeSpatial(paste0(indir,roads_file))
proj4string(roads_shp)=CRS_L93
e <- extent(344999,364998,6679996,6701004)
rc <- crop(roads_shp, e)


# Init for plot
png(paste(outdir,"Figure7a.png",sep=""),width=1200, height=1200)
par(mar=c(5.1, 6.1, 4.1, 2.1))
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15), 4, 4, byrow = TRUE))

# Loop over hours
for (i in 1:length(estim_HH_start_list)){


        estim_HH_start = estim_HH_start_list[i]
        estim_HH_end = estim_HH_end_list[i]
        estim_date3 = paste0(estim_DD,estim_MM,estim_YYYY,"_",estim_HH_start,"h",estim_HH_end,"h")
        estim_HH_start_EN = estim_HH_start_EN_list[i]
        estim_HH_end_EN = estim_HH_end_EN_list[i]
        unit <- unit_list[i]
        file_EDK = paste0("EDK_grid_Nantes_PM10_",estim_date3,"_MS_FS.csv")

# Define period
time_start = as.POSIXct(paste0(estim_date,estim_HH_start,":00:00"),origin = "1970-01-01", tz="GMT")
time_end = as.POSIXct(paste0(estim_date,estim_HH_end,":00:00"),origin = "1970-01-01", tz="GMT")

# Load estimate grid
grid_EDK <-read.csv(paste0(indir2,file_EDK),header=TRUE,sep=",",skip=0)
ras_esterr_EDK=raster(list(x=sort(unique(coordinates(grid_EDK)[,1])),y=sort(unique(coordinates(grid_EDK)[,2])),z=matrix(grid_EDK$StDev_EDK,nrow=length(sort(unique(coordinates(grid_EDK)[,1]))),byrow=F)))

min_error=minValue(ras_esterr_EDK); print(paste0("MIN ERROR=",min_error))
max_error=maxValue(ras_esterr_EDK); print(paste0("MAX ERROR=",max_error))

### PLOT

# Color and scale for raster
d=rbind(c(38,48,132),c(61,99,174),c(114,201,195),c(220,225,30),c(240,78,34),c(133,22,24))
palette=colorRampPalette(rgb(d[,1],d[,2],d[,3],maxColorValue=255))(128)
mycol=palette

if (pol=="PM10"){ zl <- c(4,8)}
ras_esterr_EDK[which(ras_esterr_EDK@data@values < zl[1])] <- zl[1]
ras_esterr_EDK[which(ras_esterr_EDK@data@values > zl[2])] <- zl[2]
ncol=100 # color and scale
brks    <- seq(zl[1], zl[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)

plot(ras_esterr_EDK, col=mycol_ncol, zlim=zl, breaks=brks, interpolate=FALSE, maxpixels=500000000,
     xaxt="n",yaxt="n", xaxs="i", yaxs="i",
     main=paste0(estim_HH_start_EN,'-',estim_HH_end_EN,unit), cex.main=2,
     legend=FALSE)
plot(rc,add=TRUE,col="grey50",alpha=0.3)

}
dev.off()

png(paste0(outdir,"Figure7b.png"),width=1200, height=1200)
plot.new()
image.plot(legend.only=TRUE, horizontal = TRUE, zlim=zl,breaks=brks,
       col=mycol_ncol, legend.shrink=1,legend.width=2.5,legend.mar=50,
       legend.args=list(text=bquote(.(pol) ~ (mu*g/m^3)),cex=1.5,line=1,font=1,col='black'),
       axis.args=list(at=brks[indbrks],labels=format(brkslab[indbrks],digits=2),cex.axis=1.5,col.axis="black",col="black"))
dev.off()




