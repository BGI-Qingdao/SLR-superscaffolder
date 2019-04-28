# stLFR Scaffold Assembler

## <a name=intro>Introduction</a>

  This is a scaffold assembler designed for stLFR reads[1].
It uses the link-reads information from stLFR reads to assemble contigs to scaffolds.

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
- [Contact](#contact)

## <a name=user-guide>User's Guide</a> 

### <a name=install>Installation</a>

- How to download the source codes.

*Clone all codes in 1 step.*
```
git clone --recursive  https://github.com/BGI-QingDao/stLFR_Scaffold_Assembler.git YOUR-DOWNLOAD-DIR
```
*Clone step by step*
```
git clone https://github.com/BGI-QingDao/stLFR_Scaffold_Assembler.git YOUR-DOWNLOAD-DIR
cd YOUR-DOWNLOAD-DIR
git submodule init
git submodule update
```
- How to compile source codes and install executable files.
```
cd YOUR-DOWNLOAD-DIR
./install.sh YOUR-INSTALL-DIR
```
*Notice: the intall.sh will create a new folder named by YOUR-INSTALL-DIR.*

- Structure of files in YOUR-INSTALL-DIR

```
├── scaffold                    # The pipeline script folder.
│   ├── __common_function.sh    # The utils script.
│   ├── clean_prepare.m4        # The template clean script.
│   ├── conf.ini                # The default configuration file.
│   ├── prepare.sh              # The preparation script that generates 1-step executable bash script according to template script and configuration file.
│   ├── run.sh                  # The template 1-step executable script.
│   ├── step_1_prepare_info.m4  # The template step 1 script. Run bwa-mem and parse the sam format mapping information into a series of files
│   ├── step_2_order.m4         # The template step 2 script. Determine the order of contigs to generate basic scaffolds.
│   ├── step_3_orientation.m4   # The template step 3 script. Determine the orientation those contigs in above scaffolds.
│   ├── step_4_pe_fill.m4       # The template step 4 script. Fill short contigs into above scaffolds.
│   ├── step_5_gapsize.m4       # The template step 5 script. Estimate gap sizes between those contigs in above scaffolds.
│   └── step_6_gen_seq.m4       # The template step 6 script. Generate the scaffold sequences based on all above information.
└── bin                         # The bin folder. Below utils are all binary executable files called by above scripts.
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
└── split_barcode               # a backup script for splitting reads and obtain barcodes
    ├── barcode_list.txt        # the barcode whitelist of stLFR reads
    ├── split_barcode.pl        # the splitting main script
    └── split_barcode.sh        # a easy-to-use swap of split_barcode.pl
```

### <a name=pre>Preliminary</a> 

Two input data are required: the stLFR reads and the contigs.

- the stLFR reads are required as 2 files :
    - your-prefix.read1.your-suffix
    - your-prefix.read2.your-suffix

  *We assume your stLFR reads and barcode information have already been splitted.*
  
  *If your data have not been splitted yet, then use the split barcode script below:*
```
# if your raw stLFR reads contain more than 1 lane, you need to cat all lines into a single file first!
YOUR-INSTALL-DIR/split_barcode/split_barcode.sh raw_read1.fq.gz raw_read2.fq.gz
```
*Also, you can try "1.fq_BarcodeSplit" step from stLFR_v1(https://github.com/MGI-tech-bioinformatics/stLFR_v1.git)*


**The 1st line of the resulting fastq file after read splitting should look like:**

```
@CL100050407L1C002R064_8855#514_1207_1392/1     7       1
```
*the "#514_1207_1392" part is the digital representation of original barcode sequence.*

- the contig must follow the SOAPdenovo contig format, which contains 2 files :
    - your-prefix.contig
    - your-prefix.ContigIndex

*Contigs assembled using various assemblers are acceptable.  We recommend MaSuRCA[2].*

*If your contigs generated by an assembler other than SOAPdenovo series, then you can easily convert the contigs to SOAPdenovo format by running:*

```
YOUR-INSTALL-DIR/bin/FakeSOAPContig < your-contig-sequence-file 1>your-prefix.contig 2>your-prefix.ContigIndex
```

### <a name=usage>General usage</a>

- 1st. Prepare the conf.ini file

```
cd YOUR-PROJECT-DIR
cp YOUR-INSTALL-DIR/Scaffold/conf.ini ./your-conf.ini
vim conf.ini # and configure it!
```

- 2nd. Generate the pipeline and work folder

```
YOUR-INSTALL-DIR/Scaffold/prepare.sh ./your-conf.ini
```
*This command will create a new work folder named by PROJECT_NAME from ./your-conf.ini, which has all pipeline scripts.*

*Make sure there is no folder using the same name in current path upon running this command.*

- 3rd. Run the pipeline

```
cd your-work-folder
./run.sh >log_pipeline 2>&1 &
```

*Notice : If any error happens in the middle and the running exits upon the last step, you can re-run run.sh and it will automatically detect and skip all completed previous steps.*

*Notice : Independent execution of each step is supperted. But to get the correct final output, you must re-run all subsequent steps.  That is, re-run step-[n+1 , 6] if you re-run step-n.*

### <a name=conf>Configuration</a>

```
###############################################################################
# Project Settings.
#
PROJECT_NAME="work_dir"      #name of work dirctory
THREADS=15                   #num of threads

###############################################################################
# Toos settings 
#
STLFR_ASSEMBLER_DIR=YOUR-INTALL-DIR         # stLFR Scaffold Assembler installation directory 
BWA_DIR=YOUR-BWA-DIR                        # BWA installation directory

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
SOAP_DIR="/home/soap_contig"                # the input contigs directory
SOAP_K=63                                   # the used-k-value of SOAPdenovo. If contigs are not generated using SOAPdenovo, then keep it as 63. 
PREFIX="chr19_soap2"                        # the prefix of you contig/ContigIndex

###############################################################################
# Control parameter settings 
#

###################################################
# for step 1 bwa mem
BWA_K=53                                    # the mapping k-value for bwa

###################################################
# for step 2 order
MST_BIN_SIZE=7000                           # the unit bin size for order.
MST_BIN_CLUSTER=0.1                         # the bin cluster threshold for order.
MST_CLUSTER=0.1                             # the order detect threshold.

###################################################
# for step 3 orientation
HT_BIN_SIZE=3500                            # the unit bin size for orientation.
HT_BIN_CLUSTER=0.1                          # the bin cluster threshold for orientation.
RANK=4                                      # the orientaion detect rank.

###################################################
# for step 4 pe_fill
# please make sure that SEED_MIN >= HT_BIN_SIZE * 2 
MAX_INSERT_SIZE=5000                        # the max allowed insert size.
PE_SEED_MIN=1000                            # the min contig size to fill.
PE_SEARCH_MAX=5000                          # the max search length before stop.
PE_MIN_JOINBARCODES=10                      # the min joint barcode count.
PE_MIN_COUNT=3                              # the min joint PE count.
PE_FILL=5                                   # the N size between PE-joint contigs.
###################################################
# for step 5 gap_size
GAP_BIN_SIZE=4000                           # the unit bin size for gap

###################################################
# for step 6 parameters for generate sequence 
#
MIN_SCONTIG=300                             # the min sequence length that allows to write into the final scaffolds.
MIN_N=10                                    # the min N size between 2 gapped contigs.
MIN_C=1                                     # the min N size between 2 overlaped contigs.
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

## <a name=contact>Contact</a>

- for algrothim details & discussion
    - please contact dengli1@genomics.cn or xumengyang@genomics.cn
- for code details & bug report 
    - please contact guolidong@genomics.cn 
    - or please creat an issue in github.
