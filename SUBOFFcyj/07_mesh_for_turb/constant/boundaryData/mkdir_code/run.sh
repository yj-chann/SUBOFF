#!/bin/bash
rm -rf *.e *.o 
make          
./mkdir.e   
rm -rf *.e *.o
./ofrun.sh