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

def get_subgrid_path(thID):
    currentPath = os.path.dirname(os.path.realpath(__file__))
    return currentPath+'/results/theory_'+str(thID)+'/subgrids/'

def get_target_path(thID):
    currentPath = os.path.dirname(os.path.realpath(__file__))
    return currentPath+'/results/theory_'+str(thID)+'/fastkernel/'

def get_subgrids(source, target, cur):
    cur.execute('SELECT id FROM '+source_dict[source]+' WHERE fktarget = \''+target+'\' ORDER BY id' )
    subgrids = cur.fetchall()
    if len(subgrids) <= 0:
        print("Cannot find any subgrids for target: " + target)
        exit(-1)
    return subgrids

def get_missing(subgrids, thID, target):
    missing_subgrids = []
    subgridPath = get_subgrid_path(thID) 
    for subgrid in subgrids:
        constituent = subgridPath+"FK_"+target+"_"+str(subgrid[0])+".subgrid.dat"
        if os.path.isfile(constituent) != True: 
            missing_subgrids.append(subgrid[0])
    return missing_subgrids

def get_constituents(subgrids, thID, target):
    subgridPath = get_subgrid_path(thID)
    mergeConstituents = ''
    for subgrid in subgrids:
        constituent = subgridPath+"FK_"+target+"_"+str(subgrid[0])+".subgrid.dat"
        mergeConstituents = mergeConstituents + constituent + ' '
    return mergeConstituents

def mergegrid(target, source, thID, cur):
    print("Processing grid " + bcolors.HEADER + target + bcolors.ENDC+ ", theory " + bcolors.HEADER+ str(thID) + bcolors.ENDC)
    # Search for subgrids
    currentPath = os.path.dirname(os.path.realpath(__file__))
    subgrids = get_subgrids(source, target, cur)

    mergeTarget = get_target_path(thID)+"FK_"+target+".dat"
    missing_subgrids  = get_missing(subgrids, thID, target)
    mergeConstituents = get_constituents(subgrids, thID, target)

    if len(missing_subgrids) == 0:
        os.system("FKmerge2 " +mergeTarget + ' ' + mergeConstituents)  
        print(target + bcolors.OKGREEN + " successfully generated!" + bcolors.ENDC)
    else:
        for missing in missing_subgrids:
            print(" -   Missing subgrid: " + str(missing))
        print(target + bcolors.FAIL + " generation failed!" + bcolors.ENDC)
