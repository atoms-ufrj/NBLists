#!/bin/bash
commands+=("./test_c_md data.inp")
commands+=("./test_c_md 4 data.inp")
#commands+=("")
#commands+=("")
#commands+=("")


for i in "${!commands[@]}"; do
  echo "==============================================================================="
  echo ${commands[$i]}
  echo "-------------------------------------------------------------------------------"
  eval ${commands[$i]}
  if [ "$?" == "0" ]; then
    echo -e "\n> Passed"
  else
    echo -e "\n> FAILED!"
  fi 
  echo "==============================================================================="
  echo
done

