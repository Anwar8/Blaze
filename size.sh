#!/bin/bash

# Find all source and header files and count the lines

find . -name '*.cpp' -o -name '*.hpp'  | xargs wc -l
