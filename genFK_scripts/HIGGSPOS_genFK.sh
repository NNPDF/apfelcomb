#!/bin/bash

if [ -e "FK_HIGGSPOS.dat" ] 
then
	exit 0;
else
	exit 1;
fi