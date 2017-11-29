#!/bin/bash
# This script generates the binary database from the dump
cd $(dirname $0)
if [ ! -f apfelcomb.db ]; then
    sqlite3 apfelcomb.db < apfelcomb.dat
fi
