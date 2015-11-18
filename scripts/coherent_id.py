#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 17:23:30 2015

@author: zah
"""

import argparse

import sqlite3



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Change 'id' to be sequential 
     numbers starting from zero""")
    parser.add_argument("db")
    parser.add_argument('tablename')   
    args = parser.parse_args()
    conn = sqlite3.connect(args.db)
    with conn:
        cursor = conn.cursor()
        items = cursor.execute("SELECT id FROM sets ORDER BY id").fetchall()
        for good,i in enumerate(items):
            print("Old id %d becomes %d"  % (i[0], good))
            cursor.execute("UPDATE %s SET id=%d WHERE id=%d"%(args.tablename, 
                                                              good, i[0]))
        conn.commit()
