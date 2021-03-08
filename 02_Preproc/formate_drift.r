# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/02_Preproc/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(sp) 
library(maptools)
library(gstat)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos)
library(RColorBrewer)
library(ggmap)

###########################################################################
#                   FORMATE DRIFT FOR KRIGING                             #      
#        from ADMS-Urban outputs interpolated over the domain             #
#           Created 17/09/2018 updated 15/05/2020                         # 
#        Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr       #
###########################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

print("INITIALIZATION")

city="Nantes"
pol="PM10"
unit="(ug/m3)"
indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/INPUTS/"
dirout <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/02_Preproc/figs/"
dirout2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/OUTPUTS/"
drift_file="TIN_Nantes_PM10_MOYcorr.asc"
xmin=326361.53; # x min in Lambert 93, depend on the calculation domain
ymin=6675394.0; # y min in Lambert 93, depend on the caluclation domain
res=7 # grid resolution in m
ncols= 6619; nrows= 4466 # grid dimensions
CRS_L93 <- "+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs" # Lambert 93

#####################################
#        GRID CONSTRUCTION          #
#####################################

print("GRID")

xgrid<-c() # xgrid vector
for (i in 1:ncols){if (i==1){xgrid[i]=xmin}else{xgrid[i]=xgrid[i-1]+res}}
ygrid<-c() # ygrid vector
for (i in 1:nrows){if (i==1){ygrid[i]=ymin}else{ygrid[i]=ygrid[i-1]+res}}
grid=expand.grid(xgrid,ygrid) # create grid

#####################################
#        READ MODEL OUTPUTS         #
#####################################

ras_data = raster(paste0(indir,drift_file)) # make a raster from the initial file
ras_data[ras_data==-9999] <- NA
mat_data <- as.matrix(ras_data)
mat_data = mat_data[4466:1,1:6619] # !!! MUST BE UPDATED BY THE USER !!!
mat_data = t(mat_data)
vec_data = as.vector(mat_data)

#####################################
#   ALLOCATE DATA TO THE GRID       #
#####################################

print("CREATE RASTER")

coordinates(grid)=~Var1+Var2
proj4string(grid)=CRS(CRS_L93)
grid_coords = grid@coords
grid_data <- data.frame(X=seq(1:length(grid_coords[,1])),
                        Y=seq(1:length(grid_coords[,1])),
                        drift_pol=seq(1:length(grid_coords[,1])))
grid_data$X <- grid_coords[,1]
grid_data$Y <- grid_coords[,2]
grid_data$drift_pol = vec_data

# Raster and plot
ras_pol=raster(list(x=sort(unique(coordinates(grid)[,1])),y=sort(unique(coordinates(grid)[,2])),z=matrix(grid_data$drift_pol,nrow=length(sort(unique(coordinates(grid)[,1]))),byrow=F)))

if (pol=="PM10"){ zl <- c(10,30)}
ras_pol[ which(ras_pol@data@values < zl[1]) ] <- zl[1]
ras_pol[ which(ras_pol@data@values > zl[2]) ] <- zl[2]

# Define palette
d=rbind(c(38,48,132),c(61,99,174),c(114,201,195),c(220,225,30),c(240,78,34),c(133,22,24))
palette=colorRampPalette(rgb(d[,1],d[,2],d[,3],maxColorValue=255))(128)
mycol=palette
ncol=100 # color and scale
brks    <- seq(zl[1], zl[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)

e <- extent(345000,365000, 6680000,6701000) # !!! MUST BE UPDATED BY THE USER !!!
ras_pol1_sub=crop(ras_pol, e)
crs(ras_pol1_sub) <- CRS_L93

newproj <- "+proj=longlat +datum=WGS84"
ras_pol2 <- projectRaster(ras_pol1_sub, crs=newproj)

png(paste0(dirout,"Figure1.png"),width=900, height=900)
par(mar=c(5,5,4,6)) # margin bot left top right (need space for the lgd
plot(ras_pol1_sub, col=mycol_ncol, zlim=zl, breaks=brks, interpolate=FALSE, maxpixels=500000000,
        cex.lab=1.6,
        cex.axis=1.8,
        legend.width = 1.5, legend.shrink=0.75,
        legend.args=list(text="   (ug/m3)",cex=1.8,line=1,font=1),
        axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.8))
dev.off()

print("SAVE")

# Write the output to a csv file
write.csv(grid_data, file = paste0(dirout2,"Drift_",pol,"_",city,"_",res,"m.csv"), row.names = FALSE, col.names = TRUE)




