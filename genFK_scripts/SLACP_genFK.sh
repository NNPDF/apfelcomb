#!/bin/bash

if [ -e "FK_SLACP.dat" ] 
then
	exit 0;
else
	exit 1;
fi