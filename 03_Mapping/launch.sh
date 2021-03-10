#!/bin/bash
#MSUB -r mapping_EDK
#MSUB -Q long
#MSUB -c 24
#MSUB -n 2
#MSUB -q broadwell
#MSUB -T 172800
#MSUB -o mapping_urban_scale_EDK_%I.o
#MSUB -e mapping_urban_scale_EDK_%I.e
#MSUB -A inerdrc

Rlog=/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/03_Mapping/mapping_urban_scale_EDK.Rlog

module load gnu

R CMD BATCH --no-save --no-restore mapping_urban_scale_EDK.r $Rlog 

