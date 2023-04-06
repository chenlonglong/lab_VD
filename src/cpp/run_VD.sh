#!/bin/bash

test_data=./test_data/random_pattern.txt
time ./VD.out -data_path $test_data  > exe_VD.log
