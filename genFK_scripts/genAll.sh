#!/bin/bash

# Directory setup
if [ ! -d "../compound" ]; then
  mkdir ../compound
fi

# Directory setup
if [ ! -d "../fastkernel" ]; then
  mkdir ../fastkernel
fi

# Directory setup
if [ ! -d "../cfactor" ]; then
  mkdir ../cfactor
fi

for d in */ ; do
    echo -n "* Processing $d ...." 

    # First verify inventory contents
    invFail=0
	if [ ! -e $d/inventory ]
    then
		echo -e " \x1B[91mmissing inventory\x1B[0m"
		invFail=1
    else
    	while read p; do
		  if [ ! -e $d/"FK_"$p".dat" ]; then
			echo -e " \x1B[91mmissing table(s)\x1B[0m";
			invFail=1
			break
		  fi
		done < $d/inventory
    fi

    if [ $invFail -eq 0 ]; then
    	cd $d

    	# Attempt to run genFK if present
		if [ -e genFK.sh ]; then
			/bin/bash ./genFK.sh
		   	rc=$?; 
			if [[ $rc != 0 ]]; then 
				echo -e " \x1B[91mfailed\x1B[0m"
				cd ../
				continue
			fi
	    fi

	    # if successful or no genFK is present, tidy up
		echo -e " \x1B[96msucceeded\x1B[0m"
		mv *-COMPOUND.dat ../../compound/ &> /dev/null
		mv CF* ../../cfactor/ &> /dev/null
		mv ./FK* ../../fastkernel/ &> /dev/null
		cd ../
		rm -rf ./$d
    fi
 

done
