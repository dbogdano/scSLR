#!/bin/bash

##https://unix.stackexchange.com/questions/435794/create-new-concatenated-files-of-same-name-in-multiple-directories

input_dir=$1

output_dir=$2

threads=$3

mkdir new

while read -r name; do
    find . -type f -path "$input_dir/$name" -exec samtools merge -o $output_dir -@ 4 -c -f {} 
done <file.list

