#!/bin/bash

for i in {1..3}
do
perl /Applications/RuleBender-2.2.1-osx64/BioNetGen-2.3/BNG2.pl grant_sims.bngl > out$i
mv grant_sims_nf.species run$i.species
echo "done with run $i"
done
grep -A1 "Time" * > tmp
grep -A4 "Reading list of Species" out1 > tmp2
cat tmp tmp2 > things_bound.dat
rm tmp
rm tmp2
exit 0
