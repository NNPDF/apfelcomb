tail *-$1.sh.o* -n 2 -q | awk 'NR % 2 {print 60*$5+$7+$9/60}' | ./histogram.py --max=700 --min=0 -p -b20
