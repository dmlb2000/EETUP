#!/bin/bash

USERNAME=$1
PASSWORD=""
read -s -p "PASSWORD:" PASSWORD

kepler -runkar \
	-myemsl_username "$USERNAME" \
	-myemsl_password "$PASSWORD" \
	-simquery 'proposal=45796#NWChem.CML=test_apr_17' \
	-expquery 'proposal=45796#Tag=NMR_EXP_jun_19' \
	-upmetadata 'proposal=45796#Tag=nmrpipe_test_jun_19' \
	-simloc /home/dmlb2000/nmrtools/trunk/src/sim/Quad_NMR_Simulator \
	-isoloc /home/dmlb2000/nmrtools/trunk/src/sim/oxygen17.inp.txt \
	-conv_nmrpipe_cmd '/usr/local/nmrpipe/nmrbin.linux9/bin2pipe -bo 1056 -xN 256 -xT 128 -xMODE Complex -noswap' \
	$PWD/NWChemNMR.kar
