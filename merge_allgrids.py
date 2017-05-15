#!/usr/bin/python
from mergecommon import mergegrid
import sqlite3 as lite
import sys,os
from subprocess import call

def infoSplash():
    print('Usage: ' + sys.argv[0] + " <theoryID>")
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
    missing_tables = []
    for grid in grids:
        missing = mergegrid(grid[0], grid[1], theoryID, cur)
        for table in missing:
            missing_tables.append({'source': grid[1].lower(), 'table':table})
    
    for tab in missing_tables:
        print(tab['source']+' '+str(tab['table'])) 

except lite.Error, e:
    print("Error %s:" % e.args[0])
    sys.exit(1)
finally:
    if con:
        con.close()
