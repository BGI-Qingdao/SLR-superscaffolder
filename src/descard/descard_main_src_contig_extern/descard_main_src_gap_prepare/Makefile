.PHONY: all clean bin

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11\
				-I../\
				-lz\
				-lpthread\

DEUBG_CXX  =	${CXXFLAGS} -g
RELEASE_CXX=	${CXXFLAGS}

source_cpp =	../common/files/file_reader.cpp \
		   		../biocommon/sam_bam/sam_parser.cpp\
		   		../biocommon/pair/pair_sam_parser.cpp\
				../biocommon/fasta/fasta.cpp\
				../biocommon/fastq/fastq.cpp\
				../biocommon/seq/seq.cpp\
				../biocommon/paf/PAF.cpp\
                ../biocommon/align_common/align_result.cpp\
		   		../common/files/file_writer.cpp \
		   		../common/files/gzstream.cpp \
				../common/log/log.cpp\
				../common/log/logfilter.cpp\
				../common/time/timetools.cpp\
				../common/string/stringtools.cpp\
				../stLFR/readName2Barcode.cpp\
				../stLFR/barcodeOnContig.cpp\
				../stLFR/ContigCluster.cpp\
                ../stLFR/ScaffInfo.cpp\
				../stLFR/LineGroup.cpp\
                ../stLFR/contigPairInfo.cpp\
                ../stLFR/ONT2Gap.cpp\
				../stLFR/ContigOverlap.cpp\
				../stLFR/CBB.cpp\
				../stLFR/TagId.cpp\
				../stLFR/StringIdCache.cpp\
				../common/args/argsparser.cpp\
				../soap2/contigGraph.cpp\
				../soap2/contigFasta.cpp\
				../soap2/contigType.cpp\
				../soap2/contigIndex.cpp\

source_o		= ${source_cpp:%.cpp=%.o}

.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs =Barcode2Gap Barcode2Scaffold  ReadName2Num Reads2Gap ScaffInfo2GapSeq

all :  ${jobs}




Barcode2Gap_cpp = Barcode2Gap.cpp
Barcode2Gap_o = Barcode2Gap.o
Barcode2Gap : clean ${Barcode2Gap_o} ${source_o} ../bin
	${CXX} ${Barcode2Gap_o} ${source_o} ${DEUBG_CXX} -o Barcode2Gap
	mv $@ ../bin/




Barcode2Scaffold_cpp = Barcode2Scaffold.cpp
Barcode2Scaffold_o = Barcode2Scaffold.o
Barcode2Scaffold : clean ${Barcode2Scaffold_o} ${source_o} ../bin
	${CXX} ${Barcode2Scaffold_o} ${source_o} ${DEUBG_CXX} -o Barcode2Scaffold
	mv $@ ../bin/




ReadName2Num_cpp = ReadName2Num.cpp
ReadName2Num_o = ReadName2Num.o
ReadName2Num : clean ${ReadName2Num_o} ${source_o} ../bin
	${CXX} ${ReadName2Num_o} ${source_o} ${DEUBG_CXX} -o ReadName2Num
	mv $@ ../bin/




Reads2Gap_cpp = Reads2Gap.cpp
Reads2Gap_o = Reads2Gap.o
Reads2Gap : clean ${Reads2Gap_o} ${source_o} ../bin
	${CXX} ${Reads2Gap_o} ${source_o} ${DEUBG_CXX} -o Reads2Gap
	mv $@ ../bin/




ScaffInfo2GapSeq_cpp = ScaffInfo2GapSeq.cpp
ScaffInfo2GapSeq_o = ScaffInfo2GapSeq.o
ScaffInfo2GapSeq : clean ${ScaffInfo2GapSeq_o} ${source_o} ../bin
	${CXX} ${ScaffInfo2GapSeq_o} ${source_o} ${DEUBG_CXX} -o ScaffInfo2GapSeq
	mv $@ ../bin/


jobs_o=  ${Barcode2Gap_o} ${Barcode2Scaffold_o} ${ReadName2Num_o} ${Reads2Gap_o} ${ScaffInfo2GapSeq_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}

