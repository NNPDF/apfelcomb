tail *.o* -n 2 -q | awk 'NR % 2 {print 60*$5+$7+$9/60}'
