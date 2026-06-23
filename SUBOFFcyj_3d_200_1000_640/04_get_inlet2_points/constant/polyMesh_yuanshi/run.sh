#!/bin/bash
rm -rf grid.o
rm -rf grid.e
make
./grid.e
rm -rf grid.o
rm -rf grid.e
