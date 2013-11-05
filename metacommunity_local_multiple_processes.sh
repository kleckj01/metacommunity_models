#!/bin/bash

echo "Starting 4 processes"
for i in {1..4} ;
do
octave --eval "metacommunity('diary$i.txt')" &
done


