#!/bin/bash
for i in {0..383} ; do mkdir -p "../processor$i/0"; done
for i in {0..383} ; do cp -r "../../09_3d_startup/processor$i/constant"  "../processor$i" ; done
for i in {0..383} ; do cp -r "../../09_3d_startup/processor$i/1/p"  "../processor$i/0/p"; done
for i in {0..383} ; do cp -r "../../09_3d_startup/processor$i/1/U"  "../processor$i/0/U"; done
for i in {0..383} ; do cp -r "../../09_3d_startup/processor$i/1/nut"  "../processor$i/0/nut"; done