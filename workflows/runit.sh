#!/bin/bash

USERNAME=$1
PASSWORD=""
read -s -p "PASSWORD:" PASSWORD

kepler -runkar \
	-username "$USERNAME" \
	-password "$PASSWORD" \
	-query 'proposal=45796#NWChem.CML=test_apr_17' \
	-upmetadata 'proposal=45796#Tag=nmrpipe_test_may_16' \
	-simloc /home/dmlb2000/nmrtools/trunk/src/sim/Quad_NMR_Simulator \
	-isoloc /home/dmlb2000/nmrtools/trunk/src/sim/oxygen17.inp.txt \
	$PWD/NWChemNMR.kar
