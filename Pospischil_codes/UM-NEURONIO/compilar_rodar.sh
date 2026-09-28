#!/bin/bash

for i in {20..23}; do 
  icc -fast -o x$i.x HH_RS_FS_bash.c -lm
  echo 
done

for i in {20..23}; do 
   echo $i | ./x$i.x & 
   echo 
   sleep 1
   rm x$i.x 
done

