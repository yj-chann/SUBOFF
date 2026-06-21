#!/bin/bash
rm -rf *.e *.o 
make          
./suboff_main.e   
rm -rf *.e *.o
