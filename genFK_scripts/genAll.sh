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
 
    cd $d
    if [ ! -e genFK.sh ]
    then
		echo -e " \x1B[91merror\x1B[0m"
		cd ../
    else
		/bin/bash ./genFK.sh
	   	rc=$?; 
		if [[ $rc != 0 ]]; 
		then 
			echo -e " \x1B[91mfailed\x1B[0m"
			cd ../
		else
			echo -e " \x1B[96msucceeded\x1B[0m"
			mv *-COMPOUND.dat ../../compound/ &> /dev/null
			mv CF* ../../cfactor/ &> /dev/null
			mv ./FK* ../../fastkernel/ &> /dev/null
			rm genFK.sh
			cd ../
			rm -rf ./$d
		fi
    fi

done
