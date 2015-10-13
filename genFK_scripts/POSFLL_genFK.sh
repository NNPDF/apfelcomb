#!/bin/bash

if [ -e "FK_POSFLL.dat" ] 
then
	exit 0;
else
	exit 1;
fi