#!/usr/bin/python
import sqlite3 as lite
import sys,os
from subprocess import call

def infoSplash():
    print 'Usage: ' + sys.argv[0] + " <gridID> <theoryID>"
    exit(-1)

if len(sys.argv) != 3:
    infoSplash()

gridID = sys.argv[1]
theoryID = sys.argv[2]

currentPath = os.path.dirname(os.path.realpath(__file__))
subgridPath = currentPath+'/results/theory_'+str(theoryID)+'/subgrids/'
ftargetPath = currentPath+'/results/theory_'+str(theoryID)+'/fastkernel/'
sqlite3Path = currentPath+'/db/applgrid.db'

print "Processing grid " + gridID + ", theory " + theoryID

# sqlite con
con = None
try:
    # Setup sql connection
    con = lite.connect(sqlite3Path)
    cur = con.cursor()    
    
    cur.execute('SELECT setname, name FROM grids WHERE id = \''+gridID+'\' ORDER BY id ' )
    grids = cur.fetchall()
    if len(grids) != 1:
        print "RowError: Either cannot find id: " + gridID + " or something very strange has happened"
        exit(-1)

    setname = grids[0][0] # Parent dataset
    fktname = grids[0][1] # FK merge target
    print "Processing "+ setname, fktname

    # Search for subgrids
    cur.execute('SELECT id FROM subgrids WHERE fktarget = \''+fktname+'\' ' )
    subgrids = cur.fetchall()
    if len(subgrids) <= 0:
        print "Cannot find any subgrids for target: " + fktname
        exit(-1)
    print str(len(subgrids)) + " subgrids found"

    mergeTarget = ftargetPath+"FK_"+fktname+".dat"
    mergeConstituents = ''
    for subgrid in subgrids:
        constituent = subgridPath+"FK_"+fktname+"_"+str(subgrid[0])+".subgrid.dat"
        # print constituent, os.path.isfile(constituent) 
        mergeConstituents = mergeConstituents + constituent + ' '
    os.system("FKmerge2 " + mergeConstituents+' > '+ mergeTarget)  

    # print mergeConstituents
    
except lite.Error, e:
    print "Error %s:" % e.args[0]
    sys.exit(1)
finally:
    if con:
        con.close()
