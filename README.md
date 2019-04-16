# stLFR Scaffold Assembler

## <a name=intro>Introduction</a>

  This is a scaffold assembler pipeline for stLFR reads[1].
It can use the link-reads information from stLFR reads to assembler contigs into scaffolds.

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
├── Scaffold                    # The pipeline script folder.
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
    ├── BinCluster
    ├── ChopBin
    ├── FakeSOAPContig
    ├── FillTrunkByPE
    ├── GapSize
    ├── MST
    ├── Orientation
    ├── PEGraph
    ├── ParseReadName
    ├── Sam2ReadOnContig
    ├── ScaffInfo2Seq
    ├── SeedCluster
    ├── SplitInfo
    ├── StaticsticUnique
    └── Trunk2ScaffInfo
```

### <a name=pre>Preliminary</a> 

This pipeline need two input data : the stLFR reads and the contig .

- the stLFR reads must have and only have 2 file :
    - your-prefix.read1.your-suffix
    - your-prefix.read2.your-suffix

*If you only have the raw reads of stLFR technology, before run this pipeline , you need do format convert first, see all details from :*

```
//TODO : the split barcode page
```

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
//TODO
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

