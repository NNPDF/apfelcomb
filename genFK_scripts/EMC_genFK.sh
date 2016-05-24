#!/bin/bash

if [ -e "FK_EMC.dat" ] 
then
	exit 0;
else
	exit 1;
fi
