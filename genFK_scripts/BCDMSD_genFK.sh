#!/bin/bash

if [ -e "FK_BCDMSD.dat" ] 
then
	exit 0;
else
	exit 1;
fi