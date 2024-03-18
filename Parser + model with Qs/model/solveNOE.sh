#!/bin/bash
SOLVERS="baron ipopt snopt conopt minos knitro"
for s in $SOLVERS ; do 
  echo "solving NOE using solver $s"
  ln -sf ${s}.run solverchoice.run
  time ampl NOE_ms.run
done