#!/bin/bash
# This script dumps and apfelcomb database into a text file
sqlite3 apfelcomb.db ".dump" > apfelcomb.dat
