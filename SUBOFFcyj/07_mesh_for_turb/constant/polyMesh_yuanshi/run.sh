#!/bin/bash
rm -rf *.e *.o 
make          
./grid.e   
rm -rf *.e *.o
