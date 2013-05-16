#!/bin/bash

USERNAME=$1
PASSWORD=""
read -s -p "PASSWORD:" PASSWORD

kepler -runkar \
	-username "$USERNAME" \
	-password "$PASSWORD" \
	-query 'proposal=45796#NWChem.CML=test_apr_17' \
	-upmetadata 'proposal=45796#Tag=nmrpipe_test_may_16' \
	NWChemNMR.kar
