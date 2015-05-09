#!/bin/bash

cd '/home/litao/molecule-recognition/DATA/'
if
do
nwchem NW.nw >& NW.out
then
grep -i 'scf energy' <NW.out |awk '{print $5}' >temp.txt
fi`
