#!/usr/bin/python
from __future__ import print_function
import sqlite3 as lite
import sys,os
from subprocess import call

source_dict = { 'APP': 'app_subgrids',
                'DIS': 'dis_subgrids',
                'DYP': 'dyp_subgrids'}

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def mergegrid(target, source, thID, cur):
    print("Processing grid " + bcolors.HEADER + target + bcolors.ENDC+ ", theory " + bcolors.HEADER+ str(thID) + bcolors.ENDC)
    # Search for subgrids
    currentPath = os.path.dirname(os.path.realpath(__file__))
    subgridPath = currentPath+'/results/theory_'+str(thID)+'/subgrids/'
    ftargetPath = currentPath+'/results/theory_'+str(thID)+'/fastkernel/'
    cur.execute('SELECT id FROM '+source_dict[source]+' WHERE fktarget = \''+target+'\' ORDER BY id' )
    subgrids = cur.fetchall()
    if len(subgrids) <= 0:
        print("Cannot find any subgrids for target: " + target)
        exit(-1)

    mergeTarget = ftargetPath+"FK_"+target+".dat"
    mergeConstituents = ''
    missing_subgrids = []
    for subgrid in subgrids:
        constituent = subgridPath+"FK_"+target+"_"+str(subgrid[0])+".subgrid.dat"
        if os.path.isfile(constituent) != True: 
            missing_subgrids.append(subgrid[0])
            print(" -   Missing subgrid: " + str(subgrid[0]))
        mergeConstituents = mergeConstituents + constituent + ' '

    if len(missing_subgrids) == 0:
        os.system("FKmerge2 " +mergeTarget + ' ' + mergeConstituents)  
        print(target + bcolors.OKGREEN + " successfully generated!" + bcolors.ENDC)
    else:
        print(target + bcolors.FAIL + " generation failed!" + bcolors.ENDC)
    return missing_subgrids
    # print mergeConstituents