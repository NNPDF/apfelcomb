#!/bin/bash
if [ -e "FK_TPCPR.dat" ]
then
	exit 0;
else
	exit 1;
fi
