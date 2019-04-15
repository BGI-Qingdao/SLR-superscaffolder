#!/bin/bash

if [[ ! $# -eq 1 ]] ; then 
    echo "usage : $0 path_to_install"
    echo "exit ..."
    exit 1
fi

path=$1
echo "info  :   install into $path"
bindir=$path/bin
scriptdir=$path/scripts
if [[ -e $path && ! -d $path ]] ; then 
    echo "error : $path exsit and is not a directory !!! "
    echo "exit ... "
    exit 1
fi
if [[ -e $bindir ]] ; then 
    echo "error : $bindir already exsit !!! "
    echo "exit ... "
    exit 1
fi
if [[ -e $scriptdir ]] ; then 
    echo "error : $scriptdir already exsit !!! "
    echo "exit ... "
    exit 1
fi

if [[ ! -e $path ]] ; then 
    echo "info  :   create $path now ..."
    mkdir -p $path
fi

echo "info  :   start compier ...   "
cd src/main && make >>../../log_make 2>&1 
if [[ $? != 0 ]] ; then 
    echo "error :   compier unfinish!!! please check log_make for details ."
    cd  -
    exit
fi
cd  -
cd src/tools && make >>../../log_make 2>&1  
if [[ $? != 0 ]] ; then 
    echo "error :   compier unfinish!!! please check log_make for details ."
    cd  -
    exit
fi
cd  -
echo "info  :   end compier ...   "

echo "info  :   cp executive files to $path/bin "
cp -r src/bin $path
echo "info  :   cp scripts to $path/scripts "
cp -r scripts/Scaffold $path

