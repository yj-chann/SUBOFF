#!/bin/bash
rm -rf cylinder.e
rm -rf cylinder.o
make
./cylinder.e
rm -rf cylinder.e
rm -rf cylinder.o