#!/bin/bash

USERNAME=$1
PASSWORD=""
read -s -p "PASSWORD:" PASSWORD

MY_TMPDIR=$(mktemp -d)

/home/dmlb2000/kepler-nobuild/kepler-2.4/kepler.sh -runkar \
	-myemsl_username "$USERNAME" \
	-myemsl_password "$PASSWORD" \
	-tmpdir "$MY_TMPDIR" \
	-NMRPIPE_PREFIX '/home/dmlb2000/nmr' \
	-query_server 'a3.my.emsl.pnl.gov' \
	-data_server 'a4.my.emsl.pnl.gov' \
	-simquery 'proposal=45796#NWChem.CML=test_apr_17' \
	-expquery 'proposal=45796#Tag=NMR_EXP_jun_19' \
	-upmetadata 'proposal=45796#Tag=nmrpipe_test_jun_19' \
	-simloc /home/dmlb2000/svn-repos/nmrtools/trunk/src/sim/Quad_NMR_Simulator \
	-isoloc /home/dmlb2000/svn-repos/nmrtools/trunk/src/sim/oxygen17.inp.txt \
	-conv_nmrpipe_cmd '/home/dmlb2000/nmr/nmrbin.linux9/bin2pipe -bo 1056 -xN 256 -xT 128 -xMODE Complex -noswap' \
	$PWD/NWChemNMR.kar
echo $MY_TMPDIR
ls $MY_TMPDIR
