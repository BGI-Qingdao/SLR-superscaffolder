# Table of contents in main folder

Each file in this folder will be compiled as a standalone binary executable application in bin folder.

## FakeSOAPContig

Convert other format of a input draft assembly into SOAPdenovo2's contig format.

## ParseReadName

Parse the text information of readname of each read in read1.fastq.gz to a number and statistic the frequence of barcodes as:  

* map long string read names into numbers.
* map long string barcode names into numbers.
* calulate the frequence of barcodes and mask barcodes with low frequence.

## Sam2ReadOnContig & SplitInfo

Together they parse sam file and save data in the custom format to reduce the CPU and Memory consumption in subsequenct processes.

Notice: supplemental alignments or secondary alignments will be filtered.

Notice: This two programs will be integrated in the future.

## StatisticUnique

Choose seed contigs based on the length and mapped read coverage of each contig.

## ChopBin

Chop contigs into bins with an equal size and assign barcodes into each bin.

## BinCluster

Calculate the Jaccard similarity of each pair of contigs as:

* Calculate the Jaccard similarity (JS) of each pair of bins.
* Choose the highest JS from all pairs of bins in each pair of contigs to represent the co-barcoding correlation strength of the paired contig.

## MST

Determine the order between seed contigs as:

* Construct the co-barcoding correlation graph 
* Simplify the graph by iteratively deleting the junctions in the Minimum-Spanning-Tree of the graph without tips.
* Print the branches of the Minimum-Spanning-Tree of the simplified co-barcoding correlation graph as ordered chain of contigs.

## Orientation

Determine the orientation of each contig in the already ordered chain.

## GapSize

Evaluate the gap size between two neighboring contigs linked by co-barcoding information as:

* Fit the relation between correlation strength and distance of two bins from the same contigs by a linear-regression model.
* Use the fitted relation to predict the gap size.

## SeedCluster

Cluster contigs for local scaffolding with using PE information to fill gaps.

## PEGraph && FillTrunkByPE

Construct PE scaffold graph and fill gaps by local scaffolding.

Notice: This two programs will be integrated in the future.

## Trunk2ScaffInfo && ScaffInfo2Seq

Integrate previous result and construct the final scaffolding results.

Notice: The final results include one fasta sequence file and one agp format file for recording the scaffolding information.

Notice: This two programs will be integrated in the future.

## the stLFR folder

Contain some structures special designed for this scaffolder. All of those structures are used in at least two applications.


--------------------------

Have fun
