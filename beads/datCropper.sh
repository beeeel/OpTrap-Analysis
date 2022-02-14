#!/bin/bash

mkdir originals
cd originals
mv ../*.dat .

files=`find . -name "*.dat" -print`

for f in $files
do
    echo $f
if ! [ -e $f ] 
then
    echo "File $f doesn't exist"
    exit 1
fi

if [ -e ../$f ]
then
    echo "File ../$f does exist (no overwriting!)"
    exit 1
fi

if [ $1 -gt 0 ]
then
    echo "Skipping $1 blocks from start"
    dd ibs=8000 skip=$1 if=$f of=../$f
elif [ $1 -lt 0 ]
then
   echo "Copying first $((-$1)) blocks"
   dd ibs=8000 count=$((-$1)) if=$f of=../$f
fi
done
