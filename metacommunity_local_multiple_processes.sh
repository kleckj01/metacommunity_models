#!/bin/bash

echo "Starting 4 processes"
for i in {1..4} ;
do
octave --silent --eval "metacommunity('dynamics$i.mat')" &
done


