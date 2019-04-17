# stLFR Scaffold Assembler

## <a name=intro>Introduction</a>

  This is a scaffold assemble pipeline for stLFR reads[1].
It can use the link-reads information from stLFR reads to assemble contigs into scaffolds.

## <a name=contents>Table of Contents</a>

- [Introduction](#intro)
- [Table of Contents](#contents)
- [User's Guide](#user-guide)
    - [Installation](#install)
    - [Preliminary](#pre)
    - [General usage](#usage)
    - [Configuration](#conf)
- [Miscellaneous](#misc)
- [Reference](#ref)

## <a name=user-guide>User's Guide</a> 

### <a name=install>Installation</a>

- How to download the source codes.
```
git clone https://github.com/BGI-QingDao/stLFR_Scaffold_Assembler.git YOUR-DOWNLOAD-DIR
```
- How to compiler source codes and install executable files.
```
cd YOUR-DOWNLOAD-DIR
./install.sh YOUR-INSTALL-DIR
```
*Notice: the intall.sh will create a new folder named by YOUR-INSTALL-DIR.*

- Structure of files in YOUR-INSTALL-DIR

```
├── scaffold                    # The pipeline script folder.
│   ├── __common_function.sh    # The utils script .
│   ├── clean_prepare.m4        # The template clean script .
│   ├── conf.ini                # The default configuration .
│   ├── prepare.sh              # The prepare script that generate real bash script based on template script and configuration .
│   ├── run.sh                  # The template pipeline main script .
│   ├── step_1_prepare_info.m4  # The template step 1 script. This script will run bwa and parse the sam format mapping information into a series of files
│   ├── step_2_order.m4         # The template step 2 script. Calculate the order information of contigs to generate basic scaffolds.
│   ├── step_3_orientation.m4   # The template step 3 script. Detect the orientation information of those contigs in above scaffolds.
│   ├── step_4_pe_fill.m4       # The template step 4 script. Try to fill small contigs into above scaffolds.
│   ├── step_5_gapsize.m4       # The template step 5 script. Detect the gap size between those contigs in above scaffolds.
│   └── step_6_gen_seq.m4       # The template step 6 script. Generate the scaffold sequences basic on all above information.
└── bin                         # The bin folder. Below utils are all binary executable files that called by above scripts.
│   ├── BinCluster
│   ├── ChopBin
│   ├── FakeSOAPContig
│   ├── FillTrunkByPE
│   ├── GapSize
│   ├── MST
│   ├── Orientation
│   ├── PEGraph
│   ├── ParseReadName
│   ├── Sam2ReadOnContig
│   ├── ScaffInfo2Seq
│   ├── SeedCluster
│   ├── SplitInfo
│   ├── StaticsticUnique
│   └── Trunk2ScaffInfo
└── split_barcode               # a backup script for split barcode
    ├── barcode_list.txt        # the barcode whitelist of stLFR
    ├── split_barcode.pl        # the split barcode main script
    └── split_barcode.sh        # a easy-to-use swap of split_barcode.pl
```

### <a name=pre>Preliminary</a> 

This pipeline need two input data : the stLFR reads and the contig .

- the stLFR reads must have and only have 2 file :
    - your-prefix.read1.your-suffix
    - your-prefix.read2.your-suffix

  *We assume your stLFR reads is the reads that after barcode splitted.*
*The official barcode split step is the "1.fq_BarcodeSplit" step of stLFR_v1(https://github.com/MGI-tech-bioinformatics/stLFR_v1.git)*

  *If you don't want to download the stLFR_v1 package , you can try our split barcode script:*

```
# if your raw stLFR reads contain more than 1 lane, you need to cat the into 1 file first!
YOUR-INSTALL-DIR/split_barcode/split_barcode.sh raw_read1.fq.gz raw_read2.fq.gz
```

**The 1st line of fastq reads file that after barcode splitted looks like:**

```
@CL100050407L1C002R064_8855#514_1207_1392/1     7       1
```
*the "#514_1207_1392" part is digital representation of original sequence barcode.*

- the contig must follow the SOAPdenovo contig format, which contain 2 files :
    - your-prefix.contig
    - your-prefix.ContigIndex

*Your can use any contig assembler to assemble the contig first . Our recommend contig assembler is MaSuRCA[2].*

*If you generate your contig by other contig assembler , you can easily convert your contig into SOAPdenovo format by:*

```
YOUR-INSTALL-DIR/bin/FakeSOAPContig < your-contig-sequence-file 1>your-prefix.contig 2>your-prefix.ContigIndex
```

### <a name=usage>General usage</a>

- 1st. Prepare the conf.ini

```
cd YOUR-PROJECT-DIR
cp YOUR-INSTALL-DIR/Scaffold/conf.ini ./your-conf.ini
vim conf.ini # and configure it!
```

- 2nd. Generate the pipeline and working folder

```
YOUR-INSTALL-DIR/Scaffold/prepare.sh ./your-conf.ini
```
*This command will create the new working folder that contain all pipeline scripts.*

*The working folder named by PROJECT_NAME from ./your-conf.ini*

*Make sure the working folder is not exsit before running this command.*

- 3rd. Run the pipeline

```
cd your-working-folder
./run.sh >log_pipeline 2>&1 &
```

*Notice : assume some error happend and your command exit ,you re-run this command and it will automatic skip the finished steps.*

*Notice : separate execution of each steps are supperted, but if you re-run step-n, you must also re-run step-[n+1 , 6] to get the correct output.*

### <a name=conf>Configuration</a>

```
###############################################################################
# Project Settings.
#
PROJECT_NAME="work_dir"      #the name of work dirctory.
THREADS=15                   #multi-thread-num

###############################################################################
# Toos settings 
#
STLFR_ASSEMBLER_DIR=YOUR-INTALL-DIR         # the stLFR Scaffold Assembler install directory 
BWA_DIR=YOUR-BWA-DIR                        # the bwa install directory

###############################################################################
# Input data settings
#   
#   MAKE SURE ALL PATH ARE ABSOLUTELY PATH.
#   DO NOT USE "../"  "~/" "./" !!!
#
## for input stLFR reads
R1="/home/chr19/chr19_reads1.fq.clean.gz"   # the read1 of stLFR reads
R2="/home/chr19/chr19_reads2.fq.clean.gz"   # the read2 of stLFR reads
## for input contig
SOAP_DIR="/home/soap_contig"                # the input contig directory
SOAP_K=63                                   # the used-k-value of SOAPdenovo. If you are not use SOAPdenovo to assemble the input contig , keep it as 63. 
PREFIX="chr19_soap2"                        # the prefix of you contig

###############################################################################
# Control parameter settings 
#

###################################################
# for step 1 bwa mem
BWA_K=53                                    # the mapping k-value for bwa

###################################################
# for step 2 order
MST_BIN_SIZE=7000                           # the unit bin size for order
MST_BIN_CLUSTER=0.1                         # the bin cluster threshold for order.
MST_CLUSTER=0.1                             # the order detect threshold .

###################################################
# for step 3 orientation
HT_BIN_SIZE=3500                            # the unit bin size for order
HT_BIN_CLUSTER=0.05                         # the bin cluster threshold for orientation.
RANK=4                                      # the orientaion detect rank.

###################################################
# for step 4 pe_fill
# please make sure that SEED_MIN >= HT_BIN_SIZE * 2 
MAX_INSERT_SIZE=5000                        # the max allowed insert size.
PE_SEED_MIN=1000                            # the min contig size to fill
PE_SEARCH_MAX=5000                          # the max search length before stop
PE_MIN_JOINBARCODES=10                      # the min join barcode count 
PE_MIN_COUNT=3                              # the min join PE count
PE_FILL=5                                   # the N size between PE join contigs
###################################################
# for step 5 gap_size
GAP_BIN_SIZE=4000                           # the unit bin size for gap

###################################################
# for step 6 .parameters for generate sequence 
#
MIN_SCONTIG=300                             # the min sequence length that allow to printed into the final scaffold
MIN_N=10                                    # the min N size between 2 contig
MIN_C=1                                     # the min N size when instead of have gap ,2 contig  actually have an overlap.
```
## <a name=misc>Miscellaneous</a>

- Requirements
    - Linux system && Bash
    - gcc ( with std11 support )
    - make && m4 ( default for almost every linux distribution)
- Dependency
    - We use bwa[3] to map stLFR reads against contigs.
        - you may replace this to any tools that can generate the correct SAM output.
- Limitations
    - We do hope a "not bad" contig, like N50 >= 15K and a "not bad" accuracy.
- Resources
    - The memory peak for human-like-genome-size is 150G.
    - Some steps support multi-thread :
        - bwa part .
        - BinCluster part .

## <a name=ref>Reference</a>

[1] [Efficient and unique co-barcoding of second-generation sequencing reads from long DNA molecules enabling cost effective and accurate sequencing, haplotyping, and de novo assembly][11]
 
[2] [Hybrid assembly of the large and highly repetitive genome of Aegilops tauschii, a progenitor of bread wheat, with the MaSuRCA mega-reads algorithm][22]

[3] [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM][33]

[11]: https://www.ncbi.nlm.nih.gov/pubmed/30940689 
[22]: https://genome.cshlp.org/content/27/5/787 
[33]: https://arxiv.org/abs/1303.3997

