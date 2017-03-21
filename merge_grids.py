#!/usr/bin/python
from mergecommon import mergegrid
import sqlite3 as lite
import sys,os
from subprocess import call

def infoSplash():
    print('Usage: ' + sys.argv[0] + " <gridID> <theoryID>")
    exit(-1)

if len(sys.argv) != 3:
    infoSplash()
gridID = sys.argv[1]
theoryID = sys.argv[2]

currentPath = os.path.dirname(os.path.realpath(__file__))
sqlite3Path = currentPath+'/db/apfelcomb.db'

# sqlite con
con = None
try:
    # Setup sql connection
    con = lite.connect(sqlite3Path)
    cur = con.cursor()    
    
    cur.execute('SELECT name, source FROM grids WHERE id = \''+gridID+'\' ORDER BY id ' )
    grids = cur.fetchall()
    if len(grids) != 1:
        print("RowError: Either cannot find id: " + gridID + " or something very strange has happened")
        exit(-1)

    # merge the grids

    mergegrid(grids[0][0], grids[0][1], theoryID, cur)
   
    
except lite.Error, e:
    print("Error %s:" % e.args[0])
    sys.exit(1)
finally:
    if con:
        con.close()
