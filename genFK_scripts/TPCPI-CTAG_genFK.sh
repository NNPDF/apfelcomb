#!/bin/bash
if [ -e "FK_TPCPI-CTAG.dat" ]
then
	exit 0;
else
	exit 1;
fi
