#!/bin/bash
if [ -e "FK_TPCPI.dat" ]
then
	exit 0;
else
	exit 1;
fi
