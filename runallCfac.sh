#!/bin/bash

while read p; do
  a=( $p )
  echo ${a[0]} ${a[1]}
  ./cfac_scale ${a[0]} ${a[1]}
done <allconv.dat