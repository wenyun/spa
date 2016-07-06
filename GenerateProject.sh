#!/bin/bash

set -e

BASEDIR=$(cd $(dirname $0); pwd)
CURDIR=$(pwd)

mkdir -p Makefiles; cd Makefiles;
cmake -DCMAKE_BUILD_TYPE=Release $BASEDIR

cd $CURDIR; mkdir -p Xcode; cd Xcode
cmake -G Xcode $BASEDIR
