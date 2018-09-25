#!/usr/bin/env bash

echo "Creating free energy plots..."
./fes.py fes.dat || { echo "failed to make fes plots"; exit 1; }

echo "Creating figure..."
pdflatex fes.tex || { echo "failed to make figure"; exit 1; }
pdflatex fes.tex
pdfcrop fes.pdf fes.pdf || { echo "cropping failed"; exit 1; }
