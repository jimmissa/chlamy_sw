#!/bin/bash

file=$1

while read line; do
  arr=($line)
  grep "${arr[0]}" Cr.annotation_info.txt | awk 'BEGIN{FS=OFS=ORS="\t"}{print $2}' > temp
  if [ -s temp ]; then
  echo -ne "${arr[0]}\t"
  for (( i=1; i<${#arr[@]}-1; i++ )); do
    echo -ne "${arr[i]} "
  done
  echo -ne "${arr[-1]}\t"
        cat temp | tr -d '\n'
  else
        rm temp
  fi
  echo
done < "$file" | awk '{if (length($0)>0) print $0}'

>temp
rm temp
