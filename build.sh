#!/usr/bin/env bash

if [[ $# == 0 ]] ; then

    cd src ; make ; cd -;

    mv ./src/BarcodeOnContig ./bin/
    mv ./src/BarcodeOnBin ./bin/
    mv ./src/ClusterGap ./bin/
    mv ./src/FormatBarcodeOnRef ./bin/
    mv ./src/UnicomGraph ./bin/
    mv ./src/BinCluster ./bin/
    mv ./src/ContigTypeByRef ./bin/
    mv ./src/StaticsticUnique ./bin/
    mv ./src/LinearClusterResult ./bin/
    mv ./src/Sam2ReadInContig ./bin/
    mv ./src/test ./bin/
    mv ./src/ContigGraphType ./bin/
    mv ./src/MergeClusterResult ./bin/
else
    cd src ; make $@ ; cd -
    while [[ $# -gt 0 ]]
    do
        mv ./src/$1 ./bin/
        shift
    done
fi
