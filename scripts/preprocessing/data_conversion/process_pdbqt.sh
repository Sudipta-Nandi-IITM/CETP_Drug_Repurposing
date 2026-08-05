#!/bin/bash

input_file="$1"
output_file="$2"

sed -n "$(grep -n 'MODEL\|TORSDOF' "$input_file" | cut -d: -f 1 | awk '{if(NR%2) printf "%d,",$1+1; else printf "%dp", $1}')" "$input_file" > "$output_file"

