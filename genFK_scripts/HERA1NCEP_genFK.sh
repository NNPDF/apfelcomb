#!/bin/bash

if [ -e "FK_HERA1NCEP.dat" ] 
then
	exit 0;
else
	exit 1;
fi