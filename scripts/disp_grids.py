#!/usr/bin/python

import sqlite3 as lite
import sys,os

currentPath = os.path.dirname(os.path.realpath(__file__))
variants = ["app_subgrids", "dis_subgrids", "dyp_subgrids"]
# Attempt to find tablulate
import imp
try:
    imp.find_module('tabulate')
    found = True
except ImportError:
    found = False

# Install/import tabulate
if found == False:
    os.system("pip install tabulate --user")
from tabulate import tabulate

# sqlite con
con = None

try:
    con = lite.connect(currentPath+'/../db/apfelcomb.db')
    cur = con.cursor()    
    cur.execute('SELECT SQLITE_VERSION()')
    data = cur.fetchone()
    
    for variant in variants:
        print "********** " + variant + " **********"  
        cur.execute('SELECT * FROM '+variant+' ORDER BY id' )
        col_names = [cn[0] for cn in cur.description]
        rows = cur.fetchall()
        print tabulate(rows, headers=col_names)
    
except lite.Error, e:
    print "Error %s:" % e.args[0]
    sys.exit(1)
    
finally:
    if con:
        con.close()
