#!/usr/bin/python
from mergecommon import *
import sqlite3 as lite
import sys,os
from subprocess import call

dataPath="../nnpdf/nnpdfcpp/data/"
theoryDB = dataPath+"theory.db"
pto_dict = {1:"N",2:"NN"}

def infoSplash():
    print('Usage: ' + sys.argv[0] + " [theoryID]")
    exit(-1)

# Fetch perturbative order
def get_pto(theoryID):
    con = None
    try:
        # Setup sql connection
        con = lite.connect(theoryDB)
        cur = con.cursor()
        
        cur.execute('SELECT PTO FROM TheoryIndex WHERE id =  \'' + theoryID + ' \'' )
        return cur.fetchall()[0][0]

    except lite.Error, e:
        print("Error %s:" % e.args[0])
        sys.exit(1)
    finally:
        if con:
            con.close()

if len(sys.argv) != 2:
    infoSplash()
theoryID = sys.argv[1]
ptord = get_pto(theoryID)
os.system('mkdir ./results/theory_'+theoryID+'/cfactor/')

# Find registered C-factors
sourceDir=dataPath+pto_dict[ptord]+"LOCFAC/"
registration=sourceDir+"registeredCFactors.dat"

print("Searching " + registration)
with open(registration, "r") as register:
    for entry in register:
        entry_split = entry.split()
        if entry_split[0] == "QCD":
            cmd = './src/cfac_scale ' + str(theoryID)+' '+' '.join(entry_split[1:])
            os.system(cmd)
        else:
            cmd = 'cp '+sourceDir+"CF_"+entry_split[0]+"_"+entry_split[1]+".dat ./results/theory_"+theoryID+"/cfactor/"
            os.system(cmd)
