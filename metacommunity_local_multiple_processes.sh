#!/bin/bash

echo "Starting 2 processes"
for i in {1..2} ;
do
octave --silent --eval "metacommunity('dynamics$i.mat')" &
done


