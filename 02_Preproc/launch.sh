#!/bin/bash
#MSUB -r correlation_drift_data 
#MSUB -c 24 
#MSUB -q broadwell
#MSUB -T 43200
#MSUB -o correlation_data_model_%I.o
#MSUB -e correlation_data_model_%I.e
#MSUB -A inerdrc

Rlog=/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/02_Preproc/correlation_data_model.Rlog

R CMD BATCH --no-save --no-restore correlation_data_drift.r $Rlog 

