# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/02_Preproc/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries
library(maptools)
library(RColorBrewer)
library(fields)
library(raster)
library(rgdal)
library(gstat)
library(rgeos)
library(sp)
library(spacetime)
library(lubridate)
library(chron)

#################################################################################
#                   PREPROCESSING FOR RAW MOBILE SENSOR DATA                    #
#           1) Smoothed data if necessary => too much noise: average            #
#           2) Correct the daily variation of the background concentrations     #
#                      Created 21/08/2018 updated 14/05/2020                    #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

print("INITIALIZATION")

indir <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/INPUTS/" # path for input directory
dirout <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/02_Preproc/figs/" # path for output directory plot
dirout2 <- "/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/OUTPUTS/" # path for output directory files 
pol="PM10"
Msens_file <- "dataout_mobile_atmotrack_novembre_temp_hum_press.csv" # sensor measurements file name
ref_file <- paste0("mesures_ref_qthourly_",pol,".txt") # reference measurements file name
code_sta_pol <- c(188,140,238,107,239) # reference station's code
sta_name <- c("Bouteillerie","Victor Hugo") # reference station's name
lat1 = 47.223253; lon1 = -1.539075 # Station N°1
lat2 = 47.204264; lon2 = -1.552927 # Station N°2
lat3 = 47.194049; lon3 = -1.581454 # Station N°3
lat4 = 47.252721; lon4 = -1.574653 # Station N°4
lat5 = 47.185581; lon5 = -1.591332 # Station N°5
estim_period <- "112018"
estim_year <- substr(estim_period,3,6)
estim_month <- substr(estim_period,1,2)
day_list <- sprintf("%02d", 1:30)

#####################################
#           READ REF DATA           #
#####################################

print("READ REF DATA")

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
data_pol_ref1 <- refpol
dtparts = t(as.data.frame(strsplit(as.character(data_pol_ref1$date),' ')))
times = chron(dates=dtparts[,1],times=dtparts[,2], format=c('y-m-d','h:m:s'))
options(digits=12)
index_t = as.numeric(times) # number of days since 01/01/1970
data_pol_ref1$index_t <- index_t
data_pol_ref1=data_pol_ref1[order(data_pol_ref1$index_t),]

#####################################
#         READ SENSOR DATA          #
#####################################

print("READ SENSOR DATA")

staFile <- paste0(indir,Msens_file)
stagrid <- read.csv(file=staFile, header=TRUE, sep=";",skip=0)
stagrid = stagrid[which(stagrid$pm10 >= 4),] # Apply the repreatability criterion caluculated from 01_Analyse/repeatability.r
stagrid = stagrid[which(stagrid$pm25 >= 2),] # Apply the repreatability criterion caluculated from 01_Analyse/repeatability.r
colnames(stagrid)[which(names(stagrid) == tolower(pol))] <- "pol"
stagrid <- subset(stagrid, select=-c(temperature, humidity,pressure))
stagrid$datetime = as.POSIXct(as.character(stagrid$datetime),tz='GMT')
dtparts = t(as.data.frame(strsplit(as.character(stagrid$datetime),' ')))
times = chron(dates=dtparts[,1],times=dtparts[,2], format=c('y-m-d','h:m:s'))
options(digits=12)
index_t = as.numeric(times) # Number of days since 01/01/1970
stagrid$index_t <- index_t

#####################################
#      BACKGROUND CORRECTION        #
#####################################

print("START CORRECTION")

###  Adjustment for background concentrations: day to day differences in background concentrations
for (iday in 1:length(day_list)){
    estim_day <- day_list[iday]

# A) Estimate daily background concentration from reference stations over the sampling period
data1T1 <- data_pol_ref1
averageT1 = 1 # 1 day
timevec_out1 <- seq(from=min(unique(data1T1$index_t)/averageT1),to=max(unique(data1T1$index_t)/averageT1),by=1)*averageT1
data1T1$intervalT1 <- cut(data1T1$index_t,breaks=timevec_out1)
BG_REF <- aggregate(. ~ intervalT1,data=data1T1,mean,na.action=na.pass) # daily mean
BG_REF$TIME <- as.POSIXlt((BG_REF$index_t)*24*60*60,origin = "1970-01-01", tz="GMT")
BG_REF$TIME = as.POSIXlt(substr(as.character(BG_REF$TIME), 1, 10), tz='GMT')
BG_REF_dd <- BG_REF[which(BG_REF$TIME==as.POSIXlt(paste0(estim_year,"-",estim_month,"-",estim_day),tz='GMT')),] 

# B) Estimate backgroud concentration from each sensor observations = the 1st percentile concentration within a moving #-min window
start_date <- paste0(estim_year,"-",estim_month,"-",estim_day," 00:00:00") # Select a particular day => 21/03/2018
end_date <- paste0(estim_year,"-",estim_month,"-",estim_day," 23:59:00")
times_start = as.POSIXlt(start_date, origin="1970-01-01",tz="GMT")
times_end = as.POSIXlt(end_date, origin="1970-01-01",tz="GMT")
data2 <- stagrid[which(stagrid$datetime >= times_start & stagrid$datetime <= times_end),] # Subset dataframe
nsensd <- unique(data2$id_sensor) # Unique sensor ID

if (length(data2[,1] > 0)){ # If there is data

# Init list
nq_list <- c(0.5) # list of percentiles sensitivity test on background calculation on data

# Option smoothing data
SMOOTHING=FALSE

for (nn in 1:length(nq_list)){ # Loop over test for bgdata

    nq = nq_list[nn]

# Init output dataframes
info_preproc1 <- data.frame(id_sensor=integer(),
                  id_run=integer(),
                  deb=integer(),
                  end=integer(),
                  duration=integer(),
                  delta_tmax=integer(),
                  stringsAsFactors=FALSE)

info_preproc2 <- data.frame(id_sensor=integer(),
                  id_run=integer(),
                  deb=integer(),
                  end=integer(),
                  duration=integer(),
                  delta_tmax=integer(),
                  stringsAsFactors=FALSE)

data_preproc <- data.frame(id_sensor=integer(),
                 datetime=as.POSIXct(character()),
                 lat=numeric(),
                 lon=numeric(),
                 pol=numeric(),
                 index_t=numeric(),
                 pol_bgcorr=numeric(),
                 stringsAsFactors=FALSE)

df_preproc <- c("data_preproc1","data_preproc2","data_preproc3","data_preproc4","data_preproc5")#,"data_preproc10","data_preproc15")

for (ndf in 1:length(df_preproc)){
df_tmp <- data.frame(id_sensor=integer(),
                 datetime=as.POSIXct(character()),
                 lat=numeric(),
                 lon=numeric(),
                 pol=numeric(),
                 index_t=numeric(),
                 pol_bgcorr=numeric(),
                 stringsAsFactors=FALSE)
assign(df_preproc[ndf],df_tmp)
}

# Loop over sensor ID
for (nsens in 1:length(nsensd)){
    print(paste("SENSOR ID is ",nsensd[nsens],sep=""))

    data_sens <- data2[which(data2$id_sensor==nsensd[nsens]),] # Subset for sensor ID

    # Init variables
    run=0; count=0 ; maxdelta_t=0;
    data_run <- c(); 
    data_runf1 <- c(); data_runf2 <- c(); data_runf3 <- c(); data_runf4 <- c(); 
    data_runf5 <- c(); 
    data_sens2 <- c()
    idx_runs <- matrix(NA,60,6)
    idx_runs2 <- matrix(NA,60,6)
    NEW_RUN=FALSE

    # Loop over sensor data
    for (dd in 1:length(data_sens[,1])){
        count=count+1

        if (NEW_RUN ==TRUE){dd=dd-1} # Start a new run

        NEW_RUN=FALSE # Continue current run

        if (count==1){
            run=run+1 
            idx_runs[run,1]=nsensd[nsens]
            idx_runs[run,2]=run
            idx_runs[run,3]=dd
            data_run=data_sens[dd,]
            data_run=cbind(run,data_run)
        }else{
            time = data_sens$datetime[dd]
            delta_t <- difftime(time,data_sens$datetime[dd-1],units='secs')
        if (delta_t <= 300){ 
            data_run[count,2:7]=data_sens[dd,]
            data_run[count,1]=run
        if (delta_t >= maxdelta_t){maxdelta_t=delta_t}
            idx_runs[run,6]=maxdelta_t
        }else{
            idx_runs[run,4]=dd-1 # store bounds of the monitoring run data
            idx_runs[run,5]=difftime(data_sens$datetime[idx_runs[run,4]],data_sens$datetime[idx_runs[run,3]],units='hours')
            idx_runs[run,6]=maxdelta_t # store max delta t between two measurements over the monitoring run 

            print(paste("END RUN ",run,sep=""))

            #### Smoothing data for monitoring run: test time resolution ####
            if (SMOOTHING == TRUE & idx_runs[run,5]*60 > 1){
            print("!! DO SMOOTHING !!")
                
                #Define time resolution depending on run duration    
                if(idx_runs[run,5]*60 >1 & idx_runs[run,5]*60 <=2){res_time=c(1)}
                if(idx_runs[run,5]*60 >2 & idx_runs[run,5]*60 <=3){res_time=c(1,2)}
                if(idx_runs[run,5]*60 >3 & idx_runs[run,5]*60 <=4){res_time=c(1,2,3)}
                if(idx_runs[run,5]*60 >4 & idx_runs[run,5]*60 <=5){res_time=c(1,2,3,4)}
                if(idx_runs[run,5]*60 >5){res_time=c(1,2,3,4,5)}

                for (nres in 1:length(res_time)){
                    print(paste('SMOOTHING for ',res_time[nres],seq=""))
                    if (run==1){idx_runs2[run,1]=1}else{idx_runs2[run,1]=idx_runs2[run-1,2]+1} # store run bounds for smoothed data
                    data_run_T1 <- data_run
                    averageT1 = res_time[nres]*6.94e-4 # minute convert in days
                    timevec_out1 <- seq(from=min(unique(data_run_T1$index_t)/averageT1),to=max(unique(data_run_T1$index_t)/averageT1),by=1)*averageT1
                    data_run_T1$intervalT1 <- cut(data_run_T1$index_t,breaks=timevec_out1)
                    data_run2 <- aggregate(. ~ intervalT1,data=data_run_T1,mean,na.action=na.pass)
                    data_run2$datetime <-as.POSIXlt((data_run2$index_t)*24*60*60,origin = "1970-01-01",tz = "GMT")                   
                    len=length(data_run2[,1])
                    idx_runs2[run,2]=idx_runs2[run,1]+len-1
                    data_runf = data_run2 # if smoothing apply

                    #### Underwrite function: capture baseline trends in the time-series data (1st percentile over the moving window) ####
                    twindow <-5*60  # Search window in seconds = # minutes before and after the ith measurement 
                 
                    #Init
                    x_bg_pol  <- vector(mode="numeric", length=length(data_runf[,1]))

                    for (mm in 1:length(data_runf[,1])){
                         
                         Ci_pol <- data_runf$pol[mm] # instantaneous concentration 
      
                         # time window calculation                    
                         elapsed_time <- difftime(data_runf$datetime[mm],data_runf$datetime[1], units='secs') # elapsed time since the 1st measurement          
                         if (elapsed_time < twindow){tstart = data_runf$datetime[1]}else{tstart = data_runf$datetime[mm]-twindow}
                         remaining_time <- difftime(data_runf$datetime[length(data_runf[,1])],data_runf$datetime[mm], units='secs') # remaining time to the last measurement
                         if (remaining_time < twindow){tend = data_runf$datetime[length(data_runf[,1])]}else{tend = data_runf$datetime[mm]+twindow}
                         data_window <- data_runf[which(data_runf$datetime >= tstart & data_runf$datetime <= tend),] # select data for search window
                         pol_conc <- data_window$pol
                         options(digits=4)
                         x_bg_pol[mm] = as.numeric(quantile(pol_conc, c(.01))) # 1st percentile
                         Cf_pol <- (Ci_pol - x_bg_pol[mm]) + BG_REF_dd$pol
                         if (Cf_pol<0 | is.na(Cf_pol)){Cf_pol = NA}
                         data_runf$pol_bgcorr[mm] <- Cf_pol
                     }
           
                    # Plot
                    #png(filename=paste(dirout,"MS",nsensd[nsens],"_",estim_year,estim_month,estim_day,"_run",run,'_res',res_time[nres],"ss.png",sep=""),width=1200, height=400)     
                    #plot(data_runf$datetime,data_runf$pol,type='l',col='black',ylim=c(0,120),xaxt="n",ann=FALSE,lwd=1.5)
                    #axis(side=1,labels=F)
                    #lines(data_runf$datetime,data_runf$pol_bgcorr,col='blue',lwd=1.5)
                    #lines(data_runf$datetime,x_bg_pol,col='red',lty=2,lwd=2)
                    #abline(h=BG_REF_dd$pol, v = 0, col ="gray60",lwd=2)
                    #r <- as.POSIXct(round(range(data_runf$datetime), "mins"))
                    #axis.POSIXct(1, at = seq(r[1], r[2], by = "mins" ), format = "%H:%M")
                    #grid (NULL,NULL, lty = 6) 
                    #title(main=paste('S',nsensd[nsens],"  ",estim_year,"-",estim_month,"-",estim_day,"  run",run,' time res=',res_time[nres],"mins",sep=""),cex.main=1.5,xlab="Time",ylab=paste0(pol," (ug/m3)"),cex.lab=1.5)
                    #legend("topright",legend=c("raw data", "corr data","data BG","ref BG"), col=c("black", "blue","red","gray60"),horiz = TRUE, cex=1,lty=c(1,1,2,1),lwd=2, bty="n")
                    #dev.off()

                    # Store smoothed data depending on res time every dd loop time   
                    colrun <- data_runf$run
                    data_runf = subset(data_runf, select=-c(intervalT1,run))
                    data_runfb = cbind(colrun,data_runf)
                    data_runfb$datetime = as.numeric(data_runfb$datetime)
                    data_runfc = aggregate(data_runfb[,2:8],list(data_runfb[,1]),mean)
                    names(data_runfc)[1] <- "colrun"
                    data_runfc$datetime=as.POSIXct(data_runfc$datetime,origin = "1970-01-01",tz = "GMT")
                    data_runfb$datetime=as.POSIXct(data_runfb$datetime,origin = "1970-01-01",tz = "GMT")


                    if (length(res_time)==1){
                        data_runf1 <- rbind(data_runf1,data_runfb)
                        data_runf2 <- rbind(data_runf2,data_runfc)
                        data_runf3 <- rbind(data_runf3,data_runfc)
                        data_runf4 <- rbind(data_runf4,data_runfc)
                        data_runf5 <- rbind(data_runf5,data_runfc)
                    }
                    if (length(res_time)==2){
                        if (res_time[nres] == 1){data_runf1 = rbind(data_runf1,data_runfb)}
                        if (res_time[nres] == 2){data_runf2 = rbind(data_runf2,data_runfb)
                                                 data_runf3 = rbind(data_runf3,data_runfc)
                                                 data_runf4 = rbind(data_runf4,data_runfc)
                                                 data_runf5 = rbind(data_runf5,data_runfc)}
                    }
                    if (length(res_time)==3){
                        if (res_time[nres] == 2){data_runf2 <- rbind(data_runf2,data_runfb)}
                        if (res_time[nres] == 3){data_runf3 <- rbind(data_runf3,data_runfb)
                                                 data_runf4 <- rbind(data_runf4,data_runfc)
                                                 data_runf5 <- rbind(data_runf5,data_runfc)}
                    }
                    if (length(res_time)==4){
                        if (res_time[nres] == 1){data_runf1 <- rbind(data_runf1,data_runfb)}
                        if (res_time[nres] == 2){data_runf2 <- rbind(data_runf2,data_runfb)}
                        if (res_time[nres] == 3){data_runf3 <- rbind(data_runf3,data_runfb)}
                        if (res_time[nres] == 4){data_runf4 <- rbind(data_runf4,data_runfb)
                                                 data_runf5 <- rbind(data_runf5,data_runfc)}
                    }

                    if (length(res_time)==5){
                        if (res_time[nres] == 1){data_runf1 <- rbind(data_runf1,data_runfb)}
                        if (res_time[nres] == 2){data_runf2 <- rbind(data_runf2,data_runfb)}
                        if (res_time[nres] == 3){data_runf3 <- rbind(data_runf3,data_runfb)}
                        if (res_time[nres] == 4){data_runf4 <- rbind(data_runf4,data_runfb)}
                        if (res_time[nres] == 5){data_runf5 <- rbind(data_runf5,data_runfb)}
                    }

                    }
            }else{
                    data_runf = data_run # if no smoothing
                    #### Underwrite function: capture baseline trends in the time-series data (1st percentile over the moving window) ####
                    twindow <- 15*60  # Search window in seconds = # minutes before and after the ith measurement 

                    # Init storage vectors
                    x_bg_pol  <- vector(mode="numeric", length=length(data_runf[,1]))

                    for (mm in 1:length(data_runf[,1])){
                        Ci_pol <- data_runf$pol[mm] # instantaneous concentration    

                        # time window calculation
                        elapsed_time <- difftime(data_runf$datetime[mm],data_runf$datetime[1], units='secs') #elapsed time since the 1st measurement          
                        if (elapsed_time < twindow){tstart = data_runf$datetime[1]}else{tstart = data_runf$datetime[mm]-twindow}
                        remaining_time <- difftime(data_runf$datetime[length(data_runf[,1])],data_runf$datetime[mm], units='secs') #remaining time to the last measurement
                        if (remaining_time < twindow){tend = data_runf$datetime[length(data_runf[,1])]}else{tend = data_runf$datetime[mm]+twindow}
                        data_window <- data_runf[which(data_runf$datetime >= tstart & data_runf$datetime <= tend),] #select data for search window
                        pol_conc <- data_window$pol
                        options(digits=4)
                        x_bg_pol[mm] = as.numeric(quantile(pol_conc, c(nq)))
                        Cf_pol <- (Ci_pol - x_bg_pol[mm]) + BG_REF_dd$pol
                        if (Cf_pol<0 | is.na(Cf_pol)){Cf_pol = NA}
                        data_runf$pol_bgcorr[mm] <- Cf_pol                         
                    }

                    # Plot
                    #png(filename=paste(dirout,"MS",nsensd[nsens],"_",estim_year,estim_month,estim_day,"_run",run,"_bgdata_perc_",nq,"_",pol,".png",sep=""),width=1200, height=400)  
                    #plot(data_runf$datetime,data_runf$pol,type='l',col='black',ylim=c(0,120),xaxt="n",ann=FALSE,lwd=1.5, cex.axis=1.2)
                    #lines(data_runf$datetime,data_runf$pol_bgcorr,col='blue',lwd=1.5)
                    #lines(data_runf$datetime,x_bg_pol,col='red',lty=2,lwd=2)
                    #abline(h=BG_REF_dd$pol, v = 0, col ="gray60",lwd=2)
                    #r <- as.POSIXct(round(range(data_runf$datetime), "mins"))
                    #axis.POSIXct(1, at = seq(r[1], r[2], by = "mins" ), format = "%H:%M",cex.axis=1.2)
                    #grid (NULL,NULL, lty = 6) 
                    #title(main=paste('S',nsensd[nsens],"  ",estim_year,"-",estim_month,"-",estim_day,"  run",run,sep=""),cex.main=1.6,xlab="Time",ylab=paste0(pol," (ug/m3)"),cex.lab=1.6)
                    #legend("topright",legend=c("raw data", "corr data","data BG","ref BG"), col=c("black", "blue","red","gray60"),horiz = TRUE, cex=1.6,lty=c(1,1,2,1),lwd=2, bty="n")
                    #dev.off()

             }
     

            # Store adjusted data

            print("STORE AFTER SMOOTHING")

            if (SMOOTHING==TRUE){
                data_preproc1 <- rbind(data_preproc1,data_runf1); #print("STORE PREPPROC 1 OK");
                data_preproc2 <- rbind(data_preproc2,data_runf2); #print("STORE PREPPROC 2 OK");
                data_preproc3 <- rbind(data_preproc3,data_runf3); #print("STORE PREPPROC 3 OK");
                data_preproc4 <- rbind(data_preproc4,data_runf4); #print("STORE PREPPROC 4 OK");
                data_preproc5 <- rbind(data_preproc5,data_runf5); #print("STORE PREPPROC 5 OK");
    
            }else{
                data_preproc <- rbind(data_preproc,data_runf)
            }

            # Re-init for next monitoring run
            data_run <- c()
            data_runf1 <- c(); data_runf2 <- c(); data_runf3 <- c(); data_runf4 <- c(); 
            data_runf5 <- c(); 
            count=0
            maxdelta_t=0
            print("NEW MONITORING RUN")
            NEW_RUN=TRUE
        
        }
        }

        if (dd==length(data_sens[,1])){idx_runs[run,4]=dd}
        idx_runs[run,5]=difftime(data_sens$datetime[idx_runs[run,4]],data_sens$datetime[idx_runs[run,3]],units='hours')

    } # Loop over sensor measurements
    print(idx_runs)

    # Store info preproc
    idx_runs = as.data.frame(idx_runs)
    names(idx_runs)[1] <- "id_sensor"
    names(idx_runs)[2] <- "id_run"
    names(idx_runs)[3] <- "deb"
    names(idx_runs)[4] <- "end"
    names(idx_runs)[5] <- "duration"
    names(idx_runs)[6] <- "delta_tmax"
    info_preproc1 <- rbind(info_preproc1,idx_runs)
    info_preproc1=na.omit(info_preproc1)
    
} # Loop over sensor ID

# End of the preprocessing - Save data for kriging
if (SMOOTHING==TRUE){
    data_preproc1 = subset(data_preproc1, select=-c(index_t))
    save(data_preproc1,file=paste0(dirout2,"data_preproc_MS_res1.Rda"))
    data_preproc2 = subset(data_preproc2, select=-c(index_t))
    save(data_preproc2,file=paste0(dirout2,"data_preproc_MS_res2.Rda"))
    data_preproc3 = subset(data_preproc3, select=-c(index_t))
    save(data_preproc3,file=paste0(dirout2,"data_preproc_MS_res3.Rda"))
    data_preproc4 = subset(data_preproc4, select=-c(index_t))
    save(data_preproc4,file=paste0(dirout2,"data_preproc_MS_res4.Rda"))
    data_preproc5 = subset(data_preproc5, select=-c(index_t))
    save(data_preproc5,file=paste0(dirout2,"data_preproc_MS_res5.Rda"))
}else{
    data_preproc = subset(data_preproc, select=-c(index_t))
    save(data_preproc,file=(paste0(dirout2,"data_preproc_MS_bgdata_",nq,"_",estim_year,estim_month,estim_day,".Rda")))
}
}
} # If there is data
} # Loop over day
