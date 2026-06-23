#!/bin/bash
rm -rf cylinder.o
rm -rf cylinder.e
make
./cylinder.e
rm -rf cylinder.o
rm -rf cylinder.e

rm -rf *.out
