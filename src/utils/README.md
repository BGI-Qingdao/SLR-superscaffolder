# Table of contents in utils folder

This utils folder contains serveral standalone structures/functions/macros that used in this project and also can be used in other project without modification.

## agp

Simple codes to save scaffold information in AGP format. See details of AGP format in [here in NCBI](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).

## args

Some custom macros wrap the ```getopt_long_only``` to parse command-line parameters and automaticly construct the usage.

## collection

Simply class to calculate union, intersection and Jaccard similarity for a weighted-set.

## files

Interfaces for reading and writing files.
One interface to read/write both normal files and gzip format files.

## graph

Some simple vertex-edge graph implementation include an undirected graph and a directed graph.

Support classic algorithm :

* Disjoint-Set
* Minimum-Span-Tree

Other logic functions:

* Prune tip branches iteratively.

## incr\_array

A substitute for std::vector that based on list of const-size arrays.

Use this container only when unknown but huge amount of items need to be loaded.

## interval

A simple class for interval detection like detecting contained relationship or overlap relationship.

## linear\_fitting

Linear fitting with the least square method.

## log

A simple log class with a timer.

## multithread

A simple multi-threads implementation follow the threads-pool strategy.

## string

Functions to handle std::string like trim / split.

## unittest

A simple unittest module.

## misc

Some small and useful functions like reporting errors and halting the program.
