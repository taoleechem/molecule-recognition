#!/bin/bash
cd '/home/litao/molecule-recognition/DATA/'

 ./analyze TinkerMolecules.xyz E >& TinkerMolecules.txt
wait 
grep -i 'Total Potential Energy : ' <TinkerMolecules.txt |awk '{print $5}' > temp.txt
