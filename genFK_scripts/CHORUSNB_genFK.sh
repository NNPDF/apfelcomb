#!/bin/bash

if [ -e "FK_CHORUSNB.dat" ] 
then
	exit 0;
else
	exit 1;
fi