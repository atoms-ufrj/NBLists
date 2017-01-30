#!/bin/bash
commands+=("./test_c_md data.inp")
commands+=("./test_c_md 4 data.inp")
#commands+=("")
#commands+=("")
#commands+=("")


for i in "${!commands[@]}"; do
  echo "======================================================================"
  echo ${commands[$i]}
  echo "----------------------------------------------------------------------"
  eval ${commands[$i]}
  echo "======================================================================"
  echo
done
