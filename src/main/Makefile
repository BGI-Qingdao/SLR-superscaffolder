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
		   		../common/files/file_writer.cpp \
		   		../common/files/gzstream.cpp \
				../common/log/log.cpp\
				../common/log/logfilter.cpp\
				../common/time/timetools.cpp\
				../common/string/stringtools.cpp\
				../stLFR/barcodeId.cpp\
				../stLFR/readName2Barcode.cpp\
				../stLFR/barcodeOnContig.cpp\
				../stLFR/ContigCluster.cpp\
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

jobs = Barcode2Scaffold  BinCluster ChopBin ContigDlink ExternContigByPE FakePE2BC FCRFByTF FillContigRoad FillTrunkByPE Gap_OAO LinearCDG MergeContig MergePEInfo MinTree MST ParseContig2Read ParseRead2Contig ParseReadName PEGraph Sam2ReadOnContig SeedCluster SplitInfo SplitScaffWithReads StaticsticUnique Trunk2Scaff Trunk2ScaffInfo

all :  ${jobs}




Barcode2Scaffold_cpp = Barcode2Scaffold.cpp
Barcode2Scaffold_o = Barcode2Scaffold.o
Barcode2Scaffold : clean ${Barcode2Scaffold_o} ${source_o} ../bin
	${CXX} ${Barcode2Scaffold_o} ${source_o} ${DEUBG_CXX} -o Barcode2Scaffold
	mv $@ ../bin/




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




FillTrunkByPE_cpp = FillTrunkByPE.cpp
FillTrunkByPE_o = FillTrunkByPE.o
FillTrunkByPE : clean ${FillTrunkByPE_o} ${source_o} ../bin
	${CXX} ${FillTrunkByPE_o} ${source_o} ${DEUBG_CXX} -o FillTrunkByPE
	mv $@ ../bin/




Gap_OAO_cpp = Gap_OAO.cpp
Gap_OAO_o = Gap_OAO.o
Gap_OAO : clean ${Gap_OAO_o} ${source_o} ../bin
	${CXX} ${Gap_OAO_o} ${source_o} ${DEUBG_CXX} -o Gap_OAO
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




MinTree_cpp = MinTree.cpp
MinTree_o = MinTree.o
MinTree : clean ${MinTree_o} ${source_o} ../bin
	${CXX} ${MinTree_o} ${source_o} ${DEUBG_CXX} -o MinTree
	mv $@ ../bin/




MST_cpp = MST.cpp
MST_o = MST.o
MST : clean ${MST_o} ${source_o} ../bin
	${CXX} ${MST_o} ${source_o} ${DEUBG_CXX} -o MST
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




SplitScaffWithReads_cpp = SplitScaffWithReads.cpp
SplitScaffWithReads_o = SplitScaffWithReads.o
SplitScaffWithReads : clean ${SplitScaffWithReads_o} ${source_o} ../bin
	${CXX} ${SplitScaffWithReads_o} ${source_o} ${DEUBG_CXX} -o SplitScaffWithReads
	mv $@ ../bin/




StaticsticUnique_cpp = StaticsticUnique.cpp
StaticsticUnique_o = StaticsticUnique.o
StaticsticUnique : clean ${StaticsticUnique_o} ${source_o} ../bin
	${CXX} ${StaticsticUnique_o} ${source_o} ${DEUBG_CXX} -o StaticsticUnique
	mv $@ ../bin/




Trunk2Scaff_cpp = Trunk2Scaff.cpp
Trunk2Scaff_o = Trunk2Scaff.o
Trunk2Scaff : clean ${Trunk2Scaff_o} ${source_o} ../bin
	${CXX} ${Trunk2Scaff_o} ${source_o} ${DEUBG_CXX} -o Trunk2Scaff
	mv $@ ../bin/




Trunk2ScaffInfo_cpp = Trunk2ScaffInfo.cpp
Trunk2ScaffInfo_o = Trunk2ScaffInfo.o
Trunk2ScaffInfo : clean ${Trunk2ScaffInfo_o} ${source_o} ../bin
	${CXX} ${Trunk2ScaffInfo_o} ${source_o} ${DEUBG_CXX} -o Trunk2ScaffInfo
	mv $@ ../bin/


jobs_o=  ${Barcode2Scaffold_o} ${BinCluster_o} ${ChopBin_o} ${ContigDlink_o} ${ExternContigByPE_o} ${FakePE2BC_o} ${FCRFByTF_o} ${FillContigRoad_o} ${FillTrunkByPE_o} ${Gap_OAO_o} ${LinearCDG_o} ${MergeContig_o} ${MergePEInfo_o} ${MinTree_o} ${MST_o} ${ParseContig2Read_o} ${ParseRead2Contig_o} ${ParseReadName_o} ${PEGraph_o} ${Sam2ReadOnContig_o} ${SeedCluster_o} ${SplitInfo_o} ${SplitScaffWithReads_o} ${StaticsticUnique_o} ${Trunk2Scaff_o} ${Trunk2ScaffInfo_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}

