#!/usr/bin/python
from mergecommon import *
import sqlite3 as lite
import sys,os
from subprocess import call

def infoSplash():
    print('Usage: ' + sys.argv[0] + " [run script]")
    exit(-1)

if len(sys.argv) != 2:
    infoSplash()
theoryID = sys.argv[1]

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
        subgrids = get_subgrids(grid[1], grid[0], cur)
        missing  = get_missing(subgrids, theoryID, grid[0])
        for table in missing:
            print grid[1].lower(), table, theoryID 

except lite.Error, e:
    print("Error %s:" % e.args[0])
    sys.exit(1)
finally:
    if con:
        con.close()
