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

source_o		= ${source_cpp:%.cpp=%.o}

SplitContig_cpp 	=	SplitContig.cpp
SplitContig_o   =	${SplitContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MergePEInfo_cpp 	=	MergePEInfo.cpp
MergePEInfo_o   =	${MergePEInfo_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

SplitInfo_cpp 	=	SplitInfo.cpp
SplitInfo_o   =	${SplitInfo_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ParseReadName_cpp 	=	ParseReadName.cpp
ParseReadName_o   =	${ParseReadName_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ParseRead2Contig_cpp 	=	ParseRead2Contig.cpp
ParseRead2Contig_o   =	${ParseRead2Contig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ParseContig2Read_cpp 	=	ParseContig2Read.cpp
ParseContig2Read_o   =	${ParseContig2Read_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ChopBin_cpp 	=	ChopBin.cpp
ChopBin_o   =	${ChopBin_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

SplitScaffWithReads_cpp 	=	SplitScaffWithReads.cpp
SplitScaffWithReads_o   =	${SplitScaffWithReads_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

BinCluster_cpp =	BinCluster.cpp
BinCluster_o   =	${BinCluster_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

StaticsticUnique_cpp =	StaticsticUnique.cpp
StaticsticUnique_o   =	${StaticsticUnique_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

Sam2ReadOnContig_cpp 	= Sam2ReadOnContig.cpp
Sam2ReadOnContig_o 	 	= ${Sam2ReadOnContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigDlink_cpp 		= ContigDlink.cpp
ContigDlink_o			= ${ContigDlink_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FillContigRoad_cpp		= FillContigRoad.cpp
FillContigRoad_o 		= ${FillContigRoad_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MergeContig_cpp			= MergeContig.cpp
MergeContig_o			= ${MergeContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

LinearCDG_cpp 			= LinearCDG.cpp
LinearCDG_o				= ${LinearCDG_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MST_cpp 			= MST.cpp
MST_o				= ${MST_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MinTree_cpp 			= MinTree.cpp
MinTree_o				= ${MinTree_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

PEGraph_cpp 			= PEGraph.cpp
PEGraph_o				= ${PEGraph_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

Trunk2Scaff_cpp 			= Trunk2Scaff.cpp
Trunk2Scaff_o				= ${Trunk2Scaff_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


Trunk2ScaffGap_cpp 			= Trunk2ScaffGap.cpp
Trunk2ScaffGap_o				= ${Trunk2ScaffGap_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


Gap_OAO_cpp 			= Gap_OAO.cpp
Gap_OAO_o				= ${Gap_OAO_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ExternContigByPE_cpp 			= ExternContigByPE.cpp
ExternContigByPE_o				= ${ExternContigByPE_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

SeedCluster_cpp 			= SeedCluster.cpp
SeedCluster_o				= ${SeedCluster_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FillTrunkByPE_cpp 			= FillTrunkByPE.cpp
FillTrunkByPE_o				= ${FillTrunkByPE_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FCRFByTF_cpp 			= FCRFByTF.cpp
FCRFByTF_o				= ${FCRFByTF_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FakePE2BC_cpp 			= FakePE2BC.cpp
FakePE2BC_o				= ${FakePE2BC_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

Barcode2Scaffold_cpp 			= Barcode2Scaffold.cpp
Barcode2Scaffold_o				= ${Barcode2Scaffold_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	ChopBin\
				SplitScaffWithReads\
				ParseReadName\
				ParseRead2Contig\
				ParseContig2Read\
				SplitInfo\
				SplitContig\
				MergePEInfo\
				BinCluster \
				StaticsticUnique  \
				Sam2ReadOnContig \
				ContigDlink \
				FillContigRoad \
				MergeContig\
				LinearCDG\
				MinTree\
				MST\
				PEGraph\
				Trunk2Scaff\
				Trunk2ScaffGap\
				Gap_OAO\
				ExternContigByPE\
				FCRFByTF\
				FakePE2BC\
				Barcode2Scaffold\
				FillTrunkByPE\
				SeedCluster\

jobs_o 	   =	${ChopBin_o} \
				${SplitScaffWithReads_o}\
				${ParseReadName_o}\
				${ParseRead2Contig_o}\
				${ParseContig2Read_o}\
				${SplitInfo_o}\
				${MergePEInfo_o}\
				${BinCluster_o} \
				${StaticsticUnique_o}\
				${Sam2ReadOnContig_o}\
				${ContigDlink_o}\
				${FillContigRoad_o}\
				${MergeContig_o}\
				${LinearCDG_o} \
				${MinTree_o}\
				${MST_o}\
				${PEGraph_o}\
				${Trunk2Scaff_o}\
				${Trunk2ScaffGap_o}\
				${Gap_OAO_o}\
				${ExternContigByPE_o}\
				${FCRFByTF_o}\
				${FakePE2BC_o}\
				${Barcode2Scaffold_o}\
				${FillTrunkByPE_o}\
				${SeedCluster_o}\


all :  ${jobs} 

MergePEInfo : clean ${MergePEInfo_o} ../bin
	${CXX} ${MergePEInfo_o} ${DEUBG_CXX} -o MergePEInfo
	mv $@ ../bin/

SplitContig : clean ${SplitContig_o} ../bin
	${CXX} ${SplitContig_o} ${DEUBG_CXX}  -o SplitContig
	mv $@ ../bin/

SplitInfo : clean ${SplitInfo_o} ../bin
	${CXX}  ${SplitInfo_o} ${DEUBG_CXX} -o SplitInfo
	mv $@ ../bin/

ChopBin : clean ${ChopBin_o} ../bin
	${CXX}  ${ChopBin_o} ${DEUBG_CXX} -o ChopBin
	mv $@ ../bin/

SplitScaffWithReads : clean ${SplitScaffWithReads_o} ../bin
	${CXX}  ${SplitScaffWithReads_o} ${DEUBG_CXX} -o SplitScaffWithReads
	mv $@ ../bin/

ParseContig2Read : clean ${ParseContig2Read_o} ../bin
	${CXX} ${ParseContig2Read_o} ${DEUBG_CXX} -o ParseContig2Read
	mv $@ ../bin/

ParseRead2Contig : clean ${ParseRead2Contig_o}  ../bin
	${CXX} ${DEUBG_CXX} ${ParseRead2Contig_o} -o ParseRead2Contig
	mv $@ ../bin/

ParseReadName : clean ${ParseReadName_o} ../bin
	${CXX}  ${ParseReadName_o} ${DEUBG_CXX} -o ParseReadName
	mv $@ ../bin/

BinCluster : clean ${BinCluster_o} ../bin
	${CXX} ${BinCluster_o} ${DEUBG_CXX}  -o BinCluster
	mv $@ ../bin/


StaticsticUnique : clean ${StaticsticUnique_o}  ../bin
	${CXX} ${StaticsticUnique_o} ${DEUBG_CXX} -o StaticsticUnique
	mv $@ ../bin/


Sam2ReadOnContig: clean ${Sam2ReadOnContig_o}  ../bin
	${CXX}  ${Sam2ReadOnContig_o} ${DEUBG_CXX} -o Sam2ReadOnContig 
	mv $@ ../bin/

ContigDlink: clean ${ContigDlink_o} ../bin
	${CXX}  ${ContigDlink_o} ${DEUBG_CXX} -o ContigDlink
	mv $@ ../bin/

FillContigRoad: clean ${FillContigRoad_o}  ../bin
	${CXX} ${FillContigRoad_o} ${DEUBG_CXX} -o FillContigRoad 
	mv $@ ../bin/

MergeContig : clean ${MergeContig_o} ../bin
	${CXX} ${MergeContig_o}  ${DEUBG_CXX} -o MergeContig
	mv $@ ../bin/

LinearCDG: clean ${LinearCDG_o} ../bin
	${CXX} ${LinearCDG_o} ${DEUBG_CXX} -o LinearCDG 
	mv $@ ../bin/

MST: clean ${MST_o} ../bin
	${CXX} ${MST_o} ${DEUBG_CXX} -o MST 
	mv $@ ../bin/

MinTree: clean ${MinTree_o} ../bin
	${CXX} ${MinTree_o} ${DEUBG_CXX} -o MinTree 
	mv $@ ../bin/

SeedCluster: clean ${SeedCluster_o} ../bin
	${CXX} ${SeedCluster_o} ${DEUBG_CXX} -o SeedCluster 
	mv $@ ../bin/

FillTrunkByPE: clean ${FillTrunkByPE_o} ../bin
	${CXX} ${FillTrunkByPE_o} ${DEUBG_CXX} -o FillTrunkByPE 
	mv $@ ../bin/

Barcode2Scaffold: clean ${Barcode2Scaffold_o} ../bin
	${CXX}  ${Barcode2Scaffold_o} ${DEUBG_CXX} -o Barcode2Scaffold 
	mv $@ ../bin/

FakePE2BC: clean ${FakePE2BC_o} ../bin
	${CXX}  ${FakePE2BC_o} ${DEUBG_CXX} -o FakePE2BC 
	mv $@ ../bin/

FCRFByTF: clean ${FCRFByTF_o} ../bin
	${CXX}  ${FCRFByTF_o} ${DEUBG_CXX} -o FCRFByTF 
	mv $@ ../bin/

ExternContigByPE: clean ${ExternContigByPE_o} ../bin
	${CXX}  ${ExternContigByPE_o} ${DEUBG_CXX} -o ExternContigByPE 
	mv $@ ../bin/

Gap_OAO: clean ${Gap_OAO_o} ../bin
	${CXX} ${Gap_OAO_o}  ${DEUBG_CXX} -o Gap_OAO 
	mv $@ ../bin/

PEGraph: clean ${PEGraph_o} ../bin
	${CXX}  ${PEGraph_o} ${DEUBG_CXX} -o PEGraph 
	mv $@ ../bin/

Trunk2Scaff: clean ${Trunk2Scaff_o} ../bin
	${CXX} ${Trunk2Scaff_o} ${DEUBG_CXX} -o Trunk2Scaff 
	mv $@ ../bin/

Trunk2ScaffGap: clean ${Trunk2ScaffGap_o} ../bin
	${CXX} ${Trunk2ScaffGap_o} ${DEUBG_CXX} -o Trunk2ScaffGap 
	mv $@ ../bin/

dirty	   =${jobs_o} ${jobs}  

../bin:
	mkdir -p ../bin

clean:
	rm -rf ${dirty}
