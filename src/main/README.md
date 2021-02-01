# Table of contents in main folder

Each file in this folder will be compile as an standalone binary executable applicatiion in bin folder.

## FakeSOAPContig

Convert other format of input genome into SOAPdenovo2's contig format.

## ParseReadName

Parse read1.fastq.gz to

* map long string read names into numbers.
* map long string barcode names into numbers.
* calulate barcode frequence and mask low quality barcode.

## Sam2ReadOnContig & SplitInfo

Together they parse sam file and save data in custom format to reduce the CPU and Memory consumption in furture applications.

notice : supplemental alignment or secondary algnment will be filtered.

This two programs will be integrated in the future.

## StatisticUnique

Choose seed contigs based on length and read coverage of each contig.

## ChopBin

Chop the co-barcoding information of contigs into bins.

## BinCluster

Calculate the Jaccard similarity of each pair of bins.

Choose the highest JS from all pairs of bins in each pair of contigs to represent the correlation strength of them.

## MST

Construct the contig-correlation graph and simplify it.

Print the branches of Minimum-Spanning-Tree of simplified-contig-correlation graph as ordered chain of contigs.

## Orientation

Detect the orientation of each pair of contigs in the already ordered chain.

## GapSize

Fitting a linear-regression model by correlation strength of bins from same contigs.

Use the fitted model to predict the gap size in final scaffolding genome.

## SeedCluster

Cluster seed contig for local scaffolding when use PE information to fill gaps.

## PEGraph && FillTrunkByPE

Construct PE graph and fill gap by local scaffolding.

This two programs will be integrated in the future.

## Trunk2ScaffInfo && ScaffInfo2Seq

Integrate previous result and construct the final scaffolding results.

The final result include one agp format description file and one fasta sequence file.

This two programs will be integrated in the future.

## the stLFR folder

Contain some structures special designed for this scaffolder. All of thoes structures are used in at least two applications.


--------------------------

Have fun
