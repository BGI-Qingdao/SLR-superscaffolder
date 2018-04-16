.PHONY: all clean

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11 -I./\
				-I./biocommon/sam_bam\
				-I./biocommon/pair\
				-I./common/files\
				-I./common/log\
				-I./common/time\
				-I./common/args\
				-I./common/test\
				-I./common/string\
				-I./common/multithread\
				-I./soap2\
				-I./stLFR\
				-lz\
				-lpthread\

DEUBG_CXX  =	${CXXFLAGS} -g
RELEASE_CXX=	${CXXFLAGS}

source_h   =	contig_barcode.h  \
		   		common/files/file_reader.h \
		   		common/files/file_writer.h  \
		   		common/files/gzstream.h \
		   		biocommon/sam_bam/sam_parser.h\
				common/log/log.h\
				common/log/logfilter.h\
				common/time/timetools.h\
				common/test/Test.h\
				common/test/Check.h\
				common/string/stringtools.h\
				common/args/argsparser.h\

source_cpp =	contig_barcode.cpp \
		   		common/files/file_reader.cpp \
		   		common/files/file_writer.cpp \
		   		common/files/gzstream.cpp \
		   		biocommon/sam_bam/sam_parser.cpp\
		   		biocommon/pair/pair_sam_parser.cpp\
				common/log/log.cpp\
				common/log/logfilter.cpp\
				common/time/timetools.cpp\
				common/string/stringtools.cpp\
				stLFR/barcodeId.cpp\
				stLFR/readName2Barcode.cpp\
				stLFR/barcodeOnContig.cpp\
				stLFR/ContigCluster.cpp\
				stLFR/LineGroup.cpp\
				common/args/argsparser.cpp\
				soap2/loadGraph.cpp\
				soap2/contigGraph.cpp\
				soap2/contigFasta.cpp\
				soap2/contigType.cpp\



test_cpp   =	test.cpp\
				common/time/timetools_test.cpp\
				common/string/stringtools_test.cpp\
				biocommon/sam_bam/sam_parser_test.cpp\
				common/args/argsparser_test.cpp\

source_o					= ${source_cpp:%.cpp=%.o}

BarcodeOnContig_NoRef_cpp =	BarcodeOnContig_NoRef.cpp
BarcodeOnContig_NoRef_o   =	${BarcodeOnContig_NoRef_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

BarcodeOnContig_cpp =	BarcodeOnContig.cpp
BarcodeOnContig_o   =	${BarcodeOnContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

BarcodeOnBin_cpp =	BarcodeOnBin.cpp
BarcodeOnBin_o   =	${BarcodeOnBin_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FormatBarcodeOnRef_cpp =	FormatBarcodeOnRef.cpp
FormatBarcodeOnRef_o   =	${FormatBarcodeOnRef_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

BinCluster_cpp =	BinCluster.cpp
BinCluster_o   =	${BinCluster_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigTypeByRef_cpp =	ContigTypeByRef.cpp
ContigTypeByRef_o   =	${ContigTypeByRef_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigGraphType_cpp =	ContigGraphType.cpp
ContigGraphType_o   =	${ContigGraphType_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

StaticsticUnique_cpp =	StaticsticUnique.cpp
StaticsticUnique_o   =	${StaticsticUnique_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

UnicomGraph_cpp =	UnicomGraph.cpp
UnicomGraph_o   =	${UnicomGraph_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ClusterGap_cpp =	ClusterGap.cpp
ClusterGap_o   =	${ClusterGap_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_12_cpp =	job12.cpp
job_12_o   =	${job_12_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

LinearClusterResult_cpp =	LinearClusterResult.cpp
LinearClusterResult_o   =	${LinearClusterResult_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


Sam2ReadInContig_cpp 	= Sam2ReadInContig.cpp
Sam2ReadInContig_o 	 	= ${Sam2ReadInContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MergeClusterResult_cpp	= MergeClusterResult.cpp
MergeClusterResult_o 	= ${MergeClusterResult_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o} 

SuperContig_cpp 		= SuperContig.cpp
SuperContig_o			= ${SuperContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

LinkCheck_cpp			= LinkCheck.cpp
LinkCheck_o				= ${LinkCheck_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


LinkCheck_1_cpp			= LinkCheck_1.cpp
LinkCheck_1_o			= ${LinkCheck_1_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


FillContigRoad_cpp		= FillContigRoad.cpp
FillContigRoad_o 		= ${FillContigRoad_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigRoadFill2Fasta_cpp= ContigRoadFill2Fasta.cpp
ContigRoadFill2Fasta_o	= ${ContigRoadFill2Fasta_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MergeContig_cpp			= MergeContig.cpp
MergeContig_o			= ${MergeContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

SolveMulti_cpp			= SolveMulti.cpp
SolveMulti_o			= ${SolveMulti_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ClusterLengthSim_cpp 	= ClusterLengthSim.cpp
ClusterLengthSim_o		= ${ClusterLengthSim_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

test_o 	   =	${test_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	BarcodeOnContig BarcodeOnContig_NoRef FormatBarcodeOnRef BarcodeOnBin BinCluster \
				ContigTypeByRef ContigGraphType StaticsticUnique  UnicomGraph \
				ClusterGap job12 LinearClusterResult Sam2ReadInContig MergeClusterResult \
				SuperContig LinkCheck LinkCheck_1 FillContigRoad ContigRoadFill2Fasta \
				MergeContig SolveMulti ClusterLengthSim
jobs_o 	   =	${BarcodeOnContig_o} ${BarcodeOnContig_NoRef_o}  ${BarcodeOnBin_o} ${FormatBarcodeOnRef_o} \
				${BinCluster_o} ${ContigTypeByRef_o}\
				${ContigGraphType_o}  ${StaticsticUnique_o}  ${UnicomGraph_o} \
				${ClusterGap_o} ${job_12_o} ${LinearClusterResult_o} ${Sam2ReadInContig_o} \
				${MergeClusterResult_o} ${SuperContig_o} ${LinkCheck_o} \
				${FillContigRoad_o} ${ContigRoadFill2Fasta_o} ${LinkCheck_1_o} \
				${MergeContig_o} ${SolveMulti_o} ${ClusterLengthSim_o}

all : test ${jobs}



BarcodeOnContig_NoRef : clean ${BarcodeOnContig_NoRef_o}
	${CXX} ${DEUBG_CXX} ${BarcodeOnContig_NoRef_o} -o BarcodeOnContig_NoRef


BarcodeOnContig : clean ${BarcodeOnContig_o}
	${CXX} ${DEUBG_CXX} ${BarcodeOnContig_o} -o BarcodeOnContig
	
FormatBarcodeOnRef : clean ${FormatBarcodeOnRef_o}
	${CXX} ${DEUBG_CXX} ${FormatBarcodeOnRef_o} -o FormatBarcodeOnRef 

BarcodeOnBin : clean ${BarcodeOnBin_o}
	${CXX} ${DEUBG_CXX} ${BarcodeOnBin_o} -o BarcodeOnBin

BinCluster : clean ${BinCluster_o}
	${CXX} ${DEUBG_CXX} ${BinCluster_o} -o BinCluster

ContigTypeByRef : clean ${ContigTypeByRef_o}
	${CXX} ${DEUBG_CXX} ${ContigTypeByRef_o} -o ContigTypeByRef

ContigGraphType : clean ${ContigGraphType_o}
	${CXX} ${DEUBG_CXX} ${ContigGraphType_o} -o ContigGraphType

StaticsticUnique : clean ${StaticsticUnique_o}
	${CXX} ${DEUBG_CXX} ${StaticsticUnique_o} -o StaticsticUnique

UnicomGraph : clean ${UnicomGraph_o}
	${CXX} ${DEUBG_CXX} ${UnicomGraph_o} -o UnicomGraph

ClusterGap : clean ${ClusterGap_o}
	${CXX} ${DEUBG_CXX} ${ClusterGap_o} -o ClusterGap

job12 : clean ${job_12_o}
	${CXX} ${DEUBG_CXX} ${job_12_o} -o job12

LinearClusterResult : clean ${LinearClusterResult_o}
	${CXX} ${DEUBG_CXX} ${LinearClusterResult_o} -o LinearClusterResult

Sam2ReadInContig: clean ${Sam2ReadInContig_o}
	${CXX} ${DEUBG_CXX} ${Sam2ReadInContig_o} -o Sam2ReadInContig 

MergeClusterResult : clean ${MergeClusterResult_o}
	${CXX} ${DEUBG_CXX} ${MergeClusterResult_o} -o MergeClusterResult

SuperContig: clean ${SuperContig_o}
	${CXX} ${DEUBG_CXX} ${SuperContig_o} -o SuperContig

LinkCheck: clean ${LinkCheck_o}
	${CXX} ${DEUBG_CXX} ${LinkCheck_o} -o LinkCheck

LinkCheck_1: clean ${LinkCheck_1_o}
	${CXX} ${DEUBG_CXX} ${LinkCheck_1_o} -o LinkCheck_1

FillContigRoad: clean ${FillContigRoad_o}
	${CXX} ${DEUBG_CXX} ${FillContigRoad_o} -o FillContigRoad 

ContigRoadFill2Fasta : clean  ${ContigRoadFill2Fasta_o}
	${CXX} ${DEUBG_CXX} ${ContigRoadFill2Fasta_o} -o ContigRoadFill2Fasta

MergeContig : clean ${MergeContig_o}
	${CXX} ${DEUBG_CXX} ${MergeContig_o} -o MergeContig

SolveMulti : clean ${SolveMulti_o}
	${CXX} ${DEUBG_CXX} ${SolveMulti_o} -o SolveMulti

ClusterLengthSim : clean ${ClusterLengthSim_o}
	${CXX} ${DEUBG_CXX} ${ClusterLengthSim_o} -o ClusterLengthSim


test: clean ${test_o} ${source_h}
	${CXX} ${DEUBG_CXX} ${test_o} -o test


dirty	   =	${test_o} ${jobs_o} ${jobs}  test

clean:
	rm -rf ${dirty}
