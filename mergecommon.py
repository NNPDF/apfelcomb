#!/usr/bin/python
import sqlite3 as lite
import sys,os
from subprocess import call

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def mergegrid(target, thID, cur):
    print("Processing grid " + bcolors.HEADER + target + bcolors.ENDC+ ", theory " + bcolors.HEADER+ str(thID) + bcolors.ENDC)
    # Search for subgrids
    currentPath = os.path.dirname(os.path.realpath(__file__))
    subgridPath = currentPath+'/results/theory_'+str(thID)+'/subgrids/'
    ftargetPath = currentPath+'/results/theory_'+str(thID)+'/fastkernel/'
    cur.execute('SELECT id FROM subgrids WHERE fktarget = \''+target+'\' ORDER BY id' )
    subgrids = cur.fetchall()
    if len(subgrids) <= 0:
        print("Cannot find any subgrids for target: " + target)
        exit(-1)

    mergeTarget = ftargetPath+"FK_"+target+".dat"
    mergeConstituents = ''
    complete_subgrids = True
    for subgrid in subgrids:
        constituent = subgridPath+"FK_"+target+"_"+str(subgrid[0])+".subgrid.dat"
        if os.path.isfile(constituent) != True: 
            print(" -   Missing subgrid: " + str(subgrid[0]))
            complete_subgrids = False
        mergeConstituents = mergeConstituents + constituent + ' '

    if complete_subgrids == True:
        os.system("FKmerge2 " + mergeConstituents+' > '+ mergeTarget)  
        print(target + bcolors.OKGREEN + " successfully generated!" + bcolors.ENDC)
    else:
        print(target + bcolors.FAIL + " generation failed!" + bcolors.ENDC)

    # print mergeConstituents