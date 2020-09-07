#!/bin/bash
if (( $# < 1 ));then
	inputfilename="data/proc/invivo1.mat"
else
	inputfilename=$1
fi

export CUDA_VISIBLE_DEVICES=0
l_rho=0.1
l_b0=0.02
tune=1
peak_dyn=17
filename=${inputfilename/"data/proc"/"recon"}
outputfilename=$filename
echo $filename
for tdyn in $peak_dyn
do
	outputfilename=${outputfilename/".mat"/"_dyn_"$tdyn".mat"}
	python3 recon_B0.py --input $inputfilename --output $outputfilename --lambda_rho $l_rho --lambda_b0 $l_b0 --tune_b0 $tune --dyn $tdyn
done