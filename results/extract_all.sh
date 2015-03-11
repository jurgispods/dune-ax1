#!/bin/bash

grep -o "T=.*, TIT=.*, IT=.*" $1 | sed -e 's/TIT=//' -e 's/IT=//' -e 's/T=//' -e 's/,//' -e 's/,//'