#!/bin/bash
if [ -e "FK_HRSPI.dat" ]
then
	exit 0;
else
	exit 1;
fi
