#!/bin/bash

apps="Barcode2Gap\
 Barcode2Scaffold \
 BarcodeDepth \
 BinCluster\
 ChopBin\
 ContigDlink\
 CutContig\
 ExternContigByPE\
 FakePE2BC\
 FCRFByTF\
 FillContigRoad\
 FillTrunkByPE\
 FormatScaffInfo\
 Gap_OAO\
 gapfillerIDY\
 GenBinSimMatrix\
 ORGC\
 LinearCDG\
 MergeContig\
 MergePEInfo\
 MinTree\
 MST\
 IDYtool\
 IDYtool_paf\
 ONTGapCorrecter\
 ONTGapFiller\
 Overlap2Gap\
 ParseContig2Read\
 ParseRead2Contig\
 ParseReadName\
 PEGraph\
 ReadName2Num\
 Reads2Gap\
 Sam2ReadOnContig\
 ScaffInfo2Seq\
 ScaffInfo2GapSeq\
 SeedCluster\
 SplitInfo\
 SplitScaffWithReads\
 StaticsticUnique\
 SplitScaffSeq\
 Trunk2Scaff\
 Trunk2ScaffGap\
 Trunk2ScaffInfo\
"

jobs_o=" "

function GenApp()
{
local AppName=$1
jobs_o="$jobs_o \${$AppName""_o}"


echo """

$AppName"_cpp 	=	"$AppName".cpp"
$AppName"_o   =	"$AppName".o"
$AppName" : clean \${"$AppName"_o} \${source_o} ../bin"
	\${CXX} \${$AppName"_o} \${source_o} \${DEUBG_CXX}  -o "$AppName
	mv \$@ ../bin/

""">>Makefile

}

echo ".PHONY: all clean bin" >Makefile
echo """
CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11\\
				-I../\\
				-lz\\
				-lpthread\\

DEUBG_CXX  =	\${CXXFLAGS} -g
RELEASE_CXX=	\${CXXFLAGS}

source_cpp =	../common/files/file_reader.cpp \\
		   		../biocommon/sam_bam/sam_parser.cpp\\
		   		../biocommon/pair/pair_sam_parser.cpp\\
				../biocommon/fasta/fasta.cpp\\
				../biocommon/fastq/fastq.cpp\\
				../biocommon/seq/seq.cpp\\
				../biocommon/paf/PAF.cpp\\
                ../biocommon/align_common/align_result.cpp\\
		   		../common/files/file_writer.cpp \\
		   		../common/files/gzstream.cpp \\
				../common/log/log.cpp\\
				../common/log/logfilter.cpp\\
				../common/time/timetools.cpp\\
				../common/string/stringtools.cpp\\
				../stLFR/barcodeId.cpp\\
				../stLFR/readName2Barcode.cpp\\
				../stLFR/barcodeOnContig.cpp\\
				../stLFR/ContigCluster.cpp\\
                ../stLFR/ScaffInfo.cpp\\
				../stLFR/LineGroup.cpp\\
                ../stLFR/contigPairInfo.cpp\\
                ../stLFR/ONT2Gap.cpp\\
				../stLFR/ContigOverlap.cpp\\
				../stLFR/CBB.cpp\\
				../stLFR/TagId.cpp\\
				../stLFR/StringIdCache.cpp\\
				../common/args/argsparser.cpp\\
				../soap2/contigGraph.cpp\\
				../soap2/contigFasta.cpp\\
				../soap2/contigType.cpp\\
				../soap2/contigIndex.cpp\\

source_o		= \${source_cpp:%.cpp=%.o}

.cpp.o:
	\${CXX} \${DEUBG_CXX} -c \$< -o \$@

jobs =$apps

all :  \${jobs}

""" >>Makefile
for x in $apps
do
    GenApp $x
done

echo "jobs_o=$jobs_o">>Makefile
echo """
dirty	   =\${jobs_o} \${jobs} \${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf \${dirty}
""" >>Makefile 
