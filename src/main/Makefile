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

MinTree_cpp 			= MinTree.cpp
MinTree_o				= ${MinTree_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

PEGraph_cpp 			= PEGraph.cpp
PEGraph_o				= ${PEGraph_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

Trunk2Scaff_cpp 			= Trunk2Scaff.cpp
Trunk2Scaff_o				= ${Trunk2Scaff_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

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


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	ChopBin\
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
				PEGraph\
				Trunk2Scaff\
				Gap_OAO\
				ExternContigByPE\
				FCRFByTF\
				FillTrunkByPE\
				SeedCluster\

jobs_o 	   =	${ChopBin_o} \
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
				${PEGraph_o}\
				${Trunk2Scaff_o}\
				${Gap_OAO_o}\
				${ExternContigByPE_o}\
				${FCRFByTF_o}\
				${FillTrunkByPE_o}\
				${SeedCluster_o}\


all : bin  ${jobs} 

bin :
	mkdir -p ../bin

MergePEInfo : clean ${MergePEInfo_o}
	${CXX} ${DEUBG_CXX} ${MergePEInfo_o} -o MergePEInfo
	mv $@ ../bin/

SplitContig : clean ${SplitContig_o}
	${CXX} ${DEUBG_CXX} ${SplitContig_o} -o SplitContig
	mv $@ ../bin/

SplitInfo : clean ${SplitInfo_o}
	${CXX} ${DEUBG_CXX} ${SplitInfo_o} -o SplitInfo
	mv $@ ../bin/

ChopBin : clean ${ChopBin_o}
	${CXX} ${DEUBG_CXX} ${ChopBin_o} -o ChopBin
	mv $@ ../bin/

ParseContig2Read : clean ${ParseContig2Read_o}
	${CXX} ${DEUBG_CXX} ${ParseContig2Read_o} -o ParseContig2Read
	mv $@ ../bin/

ParseRead2Contig : clean ${ParseRead2Contig_o}
	${CXX} ${DEUBG_CXX} ${ParseRead2Contig_o} -o ParseRead2Contig
	mv $@ ../bin/

ParseReadName : clean ${ParseReadName_o}
	${CXX} ${DEUBG_CXX} ${ParseReadName_o} -o ParseReadName
	mv $@ ../bin/

BinCluster : clean ${BinCluster_o}
	${CXX} ${DEUBG_CXX} ${BinCluster_o} -o BinCluster
	mv $@ ../bin/


StaticsticUnique : clean ${StaticsticUnique_o}
	${CXX} ${DEUBG_CXX} ${StaticsticUnique_o} -o StaticsticUnique
	mv $@ ../bin/


Sam2ReadOnContig: clean ${Sam2ReadOnContig_o}
	${CXX} ${DEUBG_CXX} ${Sam2ReadOnContig_o} -o Sam2ReadOnContig 
	mv $@ ../bin/

ContigDlink: clean ${ContigDlink_o}
	${CXX} ${DEUBG_CXX} ${ContigDlink_o} -o ContigDlink
	mv $@ ../bin/

FillContigRoad: clean ${FillContigRoad_o}
	${CXX} ${DEUBG_CXX} ${FillContigRoad_o} -o FillContigRoad 
	mv $@ ../bin/

MergeContig : clean ${MergeContig_o}
	${CXX} ${DEUBG_CXX} ${MergeContig_o} -o MergeContig
	mv $@ ../bin/

LinearCDG: clean ${LinearCDG_o}
	${CXX} ${DEUBG_CXX} ${LinearCDG_o} -o LinearCDG 
	mv $@ ../bin/

MinTree: clean ${MinTree_o}
	${CXX} ${DEUBG_CXX} ${MinTree_o} -o MinTree 
	mv $@ ../bin/

SeedCluster: clean ${SeedCluster_o}
	${CXX} ${DEUBG_CXX} ${SeedCluster_o} -o SeedCluster 
	mv $@ ../bin/


FillTrunkByPE: clean ${FillTrunkByPE_o}
	${CXX} ${DEUBG_CXX} ${FillTrunkByPE_o} -o FillTrunkByPE 
	mv $@ ../bin/

FCRFByTF: clean ${FCRFByTF_o}
	${CXX} ${DEUBG_CXX} ${FCRFByTF_o} -o FCRFByTF 
	mv $@ ../bin/

ExternContigByPE: clean ${ExternContigByPE_o}
	${CXX} ${DEUBG_CXX} ${ExternContigByPE_o} -o ExternContigByPE 
	mv $@ ../bin/

Gap_OAO: clean ${Gap_OAO_o}
	${CXX} ${DEUBG_CXX} ${Gap_OAO_o} -o Gap_OAO 
	mv $@ ../bin/

PEGraph: clean ${PEGraph_o}
	${CXX} ${DEUBG_CXX} ${PEGraph_o} -o PEGraph 
	mv $@ ../bin/

Trunk2Scaff: clean ${Trunk2Scaff_o}
	${CXX} ${DEUBG_CXX} ${Trunk2Scaff_o} -o Trunk2Scaff 
	mv $@ ../bin/

dirty	   =${jobs_o} ${jobs}  

clean:
	rm -rf ${dirty}
