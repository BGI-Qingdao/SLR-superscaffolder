# Table of contents in utils folder

This utils folder contain serveral standalone structures/functions/macros that used in this project and also can be used in other project without modification.

## agp

Simple codes to save scaffold information in AGP format. See details of AGP format in [here in NCBI](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/).

## args

Some custom macros  swap the ```getopt_long_only``` to parsing command-line parameters and automaticly construct usage.

## collection

Simply class to calculate union, intersection and jaccard similarity for weighted-set.

## files

Interfaces for read and write file.
One interface to read/write both normal file and gzip format file.

## graph

A simple vertex-edge graph implement.

support classic algorithm :

* Depth-First-Search 
* Disjoint-Set
* Minimum-Span-Tree

other logic functions:
* iterator prune tip branches.

## incr\_array

A substitute for std::vector that based on list of const-size arrays.
Use this container only when unknown but huge amount of items need to be loaded.

## interval

A simple class for interval detection like detecting contained relationship or overlap relationship.

## linear\_fitting

Linear fitting with least sqrt method.

## log

A simple log class with a timer.

## multithread

A simple multi-threads implementation follow threads-pool strategy.

## sam

Classes for parsing SAM format. See details of SAM format by [downloading document here](https://samtools.github.io/hts-specs/SAMv1.pdf).

## seq

Classes for read/write files in FASTA/FASTQ format. I guess we all know FSATA/FASTQ format :) 

## string

Functions to handle std::string like trim / split.

## unittest

A simple unittest module.

## misc

Some small and useful functions like report error and halt program.
