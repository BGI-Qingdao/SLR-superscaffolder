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
		   		../biocommon/sam/sam_parser.cpp\
		   		../biocommon/sam/pair_sam_parser.cpp\
				../biocommon/seq/fasta.cpp\
				../biocommon/seq/fastq.cpp\
				../biocommon/seq/seq.cpp\
                ../biocommon/agp/agp.cpp\
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

jobs = BinCluster ChopBin FillTrunkByPE FakeSOAPContig GapSize MST Orientation ParseReadName PEGraph Sam2ReadOnContig ScaffInfo2Seq SeedCluster SplitInfo StaticsticUnique Trunk2ScaffInfo

all :  ${jobs}




BinCluster_cpp = BinCluster.cpp
BinCluster_o = BinCluster.o
BinCluster : clean ${BinCluster_o} ${source_o} ../bin
	${CXX} ${BinCluster_o} ${source_o} ${DEUBG_CXX} -o BinCluster
	mv $@ ../bin/




ChopBin_cpp = ChopBin.cpp
ChopBin_o = ChopBin.o
ChopBin : clean ${ChopBin_o} ${source_o} ../bin
	${CXX} ${ChopBin_o} ${source_o} ${DEUBG_CXX} -o ChopBin
	mv $@ ../bin/




FillTrunkByPE_cpp = FillTrunkByPE.cpp
FillTrunkByPE_o = FillTrunkByPE.o
FillTrunkByPE : clean ${FillTrunkByPE_o} ${source_o} ../bin
	${CXX} ${FillTrunkByPE_o} ${source_o} ${DEUBG_CXX} -o FillTrunkByPE
	mv $@ ../bin/




FakeSOAPContig_cpp = FakeSOAPContig.cpp
FakeSOAPContig_o = FakeSOAPContig.o
FakeSOAPContig : clean ${FakeSOAPContig_o} ${source_o} ../bin
	${CXX} ${FakeSOAPContig_o} ${source_o} ${DEUBG_CXX} -o FakeSOAPContig
	mv $@ ../bin/




GapSize_cpp = GapSize.cpp
GapSize_o = GapSize.o
GapSize : clean ${GapSize_o} ${source_o} ../bin
	${CXX} ${GapSize_o} ${source_o} ${DEUBG_CXX} -o GapSize
	mv $@ ../bin/




MST_cpp = MST.cpp
MST_o = MST.o
MST : clean ${MST_o} ${source_o} ../bin
	${CXX} ${MST_o} ${source_o} ${DEUBG_CXX} -o MST
	mv $@ ../bin/




Orientation_cpp = Orientation.cpp
Orientation_o = Orientation.o
Orientation : clean ${Orientation_o} ${source_o} ../bin
	${CXX} ${Orientation_o} ${source_o} ${DEUBG_CXX} -o Orientation
	mv $@ ../bin/




ParseReadName_cpp = ParseReadName.cpp
ParseReadName_o = ParseReadName.o
ParseReadName : clean ${ParseReadName_o} ${source_o} ../bin
	${CXX} ${ParseReadName_o} ${source_o} ${DEUBG_CXX} -o ParseReadName
	mv $@ ../bin/




PEGraph_cpp = PEGraph.cpp
PEGraph_o = PEGraph.o
PEGraph : clean ${PEGraph_o} ${source_o} ../bin
	${CXX} ${PEGraph_o} ${source_o} ${DEUBG_CXX} -o PEGraph
	mv $@ ../bin/




Sam2ReadOnContig_cpp = Sam2ReadOnContig.cpp
Sam2ReadOnContig_o = Sam2ReadOnContig.o
Sam2ReadOnContig : clean ${Sam2ReadOnContig_o} ${source_o} ../bin
	${CXX} ${Sam2ReadOnContig_o} ${source_o} ${DEUBG_CXX} -o Sam2ReadOnContig
	mv $@ ../bin/




ScaffInfo2Seq_cpp = ScaffInfo2Seq.cpp
ScaffInfo2Seq_o = ScaffInfo2Seq.o
ScaffInfo2Seq : clean ${ScaffInfo2Seq_o} ${source_o} ../bin
	${CXX} ${ScaffInfo2Seq_o} ${source_o} ${DEUBG_CXX} -o ScaffInfo2Seq
	mv $@ ../bin/




SeedCluster_cpp = SeedCluster.cpp
SeedCluster_o = SeedCluster.o
SeedCluster : clean ${SeedCluster_o} ${source_o} ../bin
	${CXX} ${SeedCluster_o} ${source_o} ${DEUBG_CXX} -o SeedCluster
	mv $@ ../bin/




SplitInfo_cpp = SplitInfo.cpp
SplitInfo_o = SplitInfo.o
SplitInfo : clean ${SplitInfo_o} ${source_o} ../bin
	${CXX} ${SplitInfo_o} ${source_o} ${DEUBG_CXX} -o SplitInfo
	mv $@ ../bin/




StaticsticUnique_cpp = StaticsticUnique.cpp
StaticsticUnique_o = StaticsticUnique.o
StaticsticUnique : clean ${StaticsticUnique_o} ${source_o} ../bin
	${CXX} ${StaticsticUnique_o} ${source_o} ${DEUBG_CXX} -o StaticsticUnique
	mv $@ ../bin/




Trunk2ScaffInfo_cpp = Trunk2ScaffInfo.cpp
Trunk2ScaffInfo_o = Trunk2ScaffInfo.o
Trunk2ScaffInfo : clean ${Trunk2ScaffInfo_o} ${source_o} ../bin
	${CXX} ${Trunk2ScaffInfo_o} ${source_o} ${DEUBG_CXX} -o Trunk2ScaffInfo
	mv $@ ../bin/


jobs_o=  ${BinCluster_o} ${ChopBin_o} ${FillTrunkByPE_o} ${FakeSOAPContig_o} ${GapSize_o} ${MST_o} ${Orientation_o} ${ParseReadName_o} ${PEGraph_o} ${Sam2ReadOnContig_o} ${ScaffInfo2Seq_o} ${SeedCluster_o} ${SplitInfo_o} ${StaticsticUnique_o} ${Trunk2ScaffInfo_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}

