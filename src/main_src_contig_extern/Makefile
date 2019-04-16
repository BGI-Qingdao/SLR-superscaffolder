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

jobs = ContigDlink ExternContigByPE FakePE2BC FCRFByTF FillContigRoad LinearCDG MergeContig MergePEInfo ParseContig2Read ParseRead2Contig SplitContig SplitScaffWithReads

all :  ${jobs}




ContigDlink_cpp = ContigDlink.cpp
ContigDlink_o = ContigDlink.o
ContigDlink : clean ${ContigDlink_o} ${source_o} ../bin
	${CXX} ${ContigDlink_o} ${source_o} ${DEUBG_CXX} -o ContigDlink
	mv $@ ../bin/




ExternContigByPE_cpp = ExternContigByPE.cpp
ExternContigByPE_o = ExternContigByPE.o
ExternContigByPE : clean ${ExternContigByPE_o} ${source_o} ../bin
	${CXX} ${ExternContigByPE_o} ${source_o} ${DEUBG_CXX} -o ExternContigByPE
	mv $@ ../bin/




FakePE2BC_cpp = FakePE2BC.cpp
FakePE2BC_o = FakePE2BC.o
FakePE2BC : clean ${FakePE2BC_o} ${source_o} ../bin
	${CXX} ${FakePE2BC_o} ${source_o} ${DEUBG_CXX} -o FakePE2BC
	mv $@ ../bin/




FCRFByTF_cpp = FCRFByTF.cpp
FCRFByTF_o = FCRFByTF.o
FCRFByTF : clean ${FCRFByTF_o} ${source_o} ../bin
	${CXX} ${FCRFByTF_o} ${source_o} ${DEUBG_CXX} -o FCRFByTF
	mv $@ ../bin/




FillContigRoad_cpp = FillContigRoad.cpp
FillContigRoad_o = FillContigRoad.o
FillContigRoad : clean ${FillContigRoad_o} ${source_o} ../bin
	${CXX} ${FillContigRoad_o} ${source_o} ${DEUBG_CXX} -o FillContigRoad
	mv $@ ../bin/




LinearCDG_cpp = LinearCDG.cpp
LinearCDG_o = LinearCDG.o
LinearCDG : clean ${LinearCDG_o} ${source_o} ../bin
	${CXX} ${LinearCDG_o} ${source_o} ${DEUBG_CXX} -o LinearCDG
	mv $@ ../bin/




MergeContig_cpp = MergeContig.cpp
MergeContig_o = MergeContig.o
MergeContig : clean ${MergeContig_o} ${source_o} ../bin
	${CXX} ${MergeContig_o} ${source_o} ${DEUBG_CXX} -o MergeContig
	mv $@ ../bin/




MergePEInfo_cpp = MergePEInfo.cpp
MergePEInfo_o = MergePEInfo.o
MergePEInfo : clean ${MergePEInfo_o} ${source_o} ../bin
	${CXX} ${MergePEInfo_o} ${source_o} ${DEUBG_CXX} -o MergePEInfo
	mv $@ ../bin/




ParseContig2Read_cpp = ParseContig2Read.cpp
ParseContig2Read_o = ParseContig2Read.o
ParseContig2Read : clean ${ParseContig2Read_o} ${source_o} ../bin
	${CXX} ${ParseContig2Read_o} ${source_o} ${DEUBG_CXX} -o ParseContig2Read
	mv $@ ../bin/




ParseRead2Contig_cpp = ParseRead2Contig.cpp
ParseRead2Contig_o = ParseRead2Contig.o
ParseRead2Contig : clean ${ParseRead2Contig_o} ${source_o} ../bin
	${CXX} ${ParseRead2Contig_o} ${source_o} ${DEUBG_CXX} -o ParseRead2Contig
	mv $@ ../bin/




SplitContig_cpp = SplitContig.cpp
SplitContig_o = SplitContig.o
SplitContig : clean ${SplitContig_o} ${source_o} ../bin
	${CXX} ${SplitContig_o} ${source_o} ${DEUBG_CXX} -o SplitContig
	mv $@ ../bin/




SplitScaffWithReads_cpp = SplitScaffWithReads.cpp
SplitScaffWithReads_o = SplitScaffWithReads.o
SplitScaffWithReads : clean ${SplitScaffWithReads_o} ${source_o} ../bin
	${CXX} ${SplitScaffWithReads_o} ${source_o} ${DEUBG_CXX} -o SplitScaffWithReads
	mv $@ ../bin/


jobs_o=  ${ContigDlink_o} ${ExternContigByPE_o} ${FakePE2BC_o} ${FCRFByTF_o} ${FillContigRoad_o} ${LinearCDG_o} ${MergeContig_o} ${MergePEInfo_o} ${ParseContig2Read_o} ${ParseRead2Contig_o} ${SplitContig_o} ${SplitScaffWithReads_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}

