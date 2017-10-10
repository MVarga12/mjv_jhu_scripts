#!/bin/bash

for i in {1..3}
do
./species.py run$i.species out$i.species
done
cat out*.species | sort -k1 -n > out_fin
exit 0
