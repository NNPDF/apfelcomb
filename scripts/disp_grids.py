#!/usr/bin/python

import sqlite3 as lite
import sys,os

currentPath = os.path.dirname(os.path.realpath(__file__))

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
    con = lite.connect(currentPath+'/../db/applgrid.db')
    
    cur = con.cursor()    
    cur.execute('SELECT SQLITE_VERSION()')
    
    data = cur.fetchone()
    

    print "********************** Available APPLgrids **********************"    

    cur.execute('SELECT id, setname, gridname FROM subgrids')
    col_names = [cn[0] for cn in cur.description]

    table = []

    rows = cur.fetchall()

    print tabulate(rows, headers=col_names)

    print "********************** Available CommonData **********************" 

    con = lite.connect(currentPath+'/../db/dis.db')
    cur = con.cursor() 

    cur.execute('SELECT id, setname, gridname FROM sets')
    col_names = [cn[0] for cn in cur.description]

    table = []

    rows = cur.fetchall()

    print tabulate(rows, headers=col_names)   
    
    if len(sys.argv) > 1 and sys.argv[1] == "sia":
        print "********************** Available CommonData SIA **********************" 

        con = lite.connect(currentPath+'/../db/sia.db')
        cur = con.cursor() 
        
        cur.execute('SELECT id, setname, gridname FROM sets')
        col_names = [cn[0] for cn in cur.description]

        table = []

        rows = cur.fetchall()
        
        print tabulate(rows, headers=col_names)   

    
except lite.Error, e:
    
    print "Error %s:" % e.args[0]
    sys.exit(1)
    
finally:
    
    if con:
        con.close()
