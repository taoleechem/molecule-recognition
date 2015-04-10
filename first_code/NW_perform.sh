#!/bin/bash
cd '/home/password-l/桌面/C++PRACTICE/Molecule_Recognize2/Molecule_Recognize2/DATA'
nwchem NW.nw >& NW.out
sleep 0.1
grep -i 'scf energy'<NW.out | awk '{print $5 }' >temp.txt
