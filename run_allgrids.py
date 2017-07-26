#!/usr/bin/python
from mergecommon import *
import sqlite3 as lite
import sys,os
from subprocess import call

def infoSplash():
    print('Usage: ' + sys.argv[0] + " [run script]")
    exit(-1)

if len(sys.argv) != 3:
    infoSplash()
theoryID = sys.argv[1]
runScript = sys.argv[2]

currentPath = os.path.dirname(os.path.realpath(__file__))
sqlite3Path = currentPath+'/db/apfelcomb.db'

# sqlite con
con = None
try:
    # Setup sql connection
    con = lite.connect(sqlite3Path)
    cur = con.cursor()
    
    cur.execute('SELECT name, source FROM grids ORDER BY id' )
    grids = cur.fetchall()
    for grid in grids:
        source = grid[1]
        target = grid[0]
        subgrids = get_subgrids(source, target, cur)
        missing  = get_missing(subgrids, theoryID, target)
        for table in missing:
            runcmd = runScript+' '+source.lower()+" "+str(table) + " " + str(theoryID)
   			print(target + ' ' + str(table))
            os.system(runcmd)

except lite.Error, e:
    print("Error %s:" % e.args[0])
    sys.exit(1)
finally:
    if con:
        con.close()
