.PHONY: all clean

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11\
				-I../\
				-lz\
				-lpthread\

DEUBG_CXX  =	${CXXFLAGS} -g
RELEASE_CXX=	${CXXFLAGS}

source_cpp =	../common/files/file_reader.cpp \
		   		../common/files/file_writer.cpp \
		   		../common/files/gzstream.cpp \
				../common/log/log.cpp\
				../common/log/logfilter.cpp\
				../common/time/timetools.cpp\
				../common/string/stringtools.cpp\
		   		../biocommon/sam_bam/sam_parser.cpp\
		   		../biocommon/pair/pair_sam_parser.cpp\
				../stLFR/barcodeId.cpp\
				../stLFR/readName2Barcode.cpp\
				../stLFR/barcodeOnContig.cpp\
				../stLFR/ContigCluster.cpp\
				../stLFR/LineGroup.cpp\
				../stLFR/CBB.cpp\
				../common/args/argsparser.cpp\
				../soap2/loadGraph.cpp\
				../soap2/contigGraph.cpp\
				../soap2/contigFasta.cpp\
				../soap2/contigType.cpp\
				contig_barcode.cpp


source_o		= ${source_cpp:%.cpp=%.o}

BarcodeOnContig_cpp =	BarcodeOnContig.cpp
BarcodeOnContig_o   =	${BarcodeOnContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

BarcodeOnBin_cpp =	BarcodeOnBin.cpp
BarcodeOnBin_o   =	${BarcodeOnBin_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

FormatBarcodeOnRef_cpp =	FormatBarcodeOnRef.cpp
FormatBarcodeOnRef_o   =	${FormatBarcodeOnRef_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigTypeByRef_cpp =	ContigTypeByRef.cpp
ContigTypeByRef_o   =	${ContigTypeByRef_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigGraphType_cpp =	ContigGraphType.cpp
ContigGraphType_o   =	${ContigGraphType_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

UnicomGraph_cpp =	UnicomGraph.cpp
UnicomGraph_o   =	${UnicomGraph_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ClusterGap_cpp =	ClusterGap.cpp
ClusterGap_o   =	${ClusterGap_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_12_cpp =	job12.cpp
job_12_o   =	${job_12_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

LinearClusterResult_cpp =	LinearClusterResult.cpp
LinearClusterResult_o   =	${LinearClusterResult_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

MergeClusterResult_cpp	= MergeClusterResult.cpp
MergeClusterResult_o 	= ${MergeClusterResult_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o} 

LinkCheck_cpp			= LinkCheck.cpp
LinkCheck_o				= ${LinkCheck_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

LinkCheck_1_cpp			= LinkCheck_1.cpp
LinkCheck_1_o			= ${LinkCheck_1_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ContigRoadFill2Fasta_cpp= ContigRoadFill2Fasta.cpp
ContigRoadFill2Fasta_o	= ${ContigRoadFill2Fasta_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

SolveMulti_cpp			= SolveMulti.cpp
SolveMulti_o			= ${SolveMulti_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

ClusterLengthSim_cpp 	= ClusterLengthSim.cpp
ClusterLengthSim_o		= ${ClusterLengthSim_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	BarcodeOnContig FormatBarcodeOnRef BarcodeOnBin \
				ContigTypeByRef ContigGraphType UnicomGraph \
				ClusterGap job12 LinearClusterResult MergeClusterResult \
				LinkCheck LinkCheck_1 ContigRoadFill2Fasta \
				SolveMulti ClusterLengthSim 
jobs_o 	   =	${BarcodeOnContig_o}  ${BarcodeOnBin_o} ${FormatBarcodeOnRef_o} \
				${ContigTypeByRef_o}\
				${ContigGraphType_o}   ${UnicomGraph_o} \
				${ClusterGap_o} ${job_12_o} ${LinearClusterResult_o}\
				${MergeClusterResult_o} ${LinkCheck_o} \
				${ContigRoadFill2Fasta_o} ${LinkCheck_1_o} \
				${SolveMulti_o} ${ClusterLengthSim_o} \

all : ${jobs}


BarcodeOnContig : clean ${BarcodeOnContig_o}
	${CXX} ${DEUBG_CXX} ${BarcodeOnContig_o} -o BarcodeOnContig
	
FormatBarcodeOnRef : clean ${FormatBarcodeOnRef_o}
	${CXX} ${DEUBG_CXX} ${FormatBarcodeOnRef_o} -o FormatBarcodeOnRef 

BarcodeOnBin : clean ${BarcodeOnBin_o}
	${CXX} ${DEUBG_CXX} ${BarcodeOnBin_o} -o BarcodeOnBin

ContigTypeByRef : clean ${ContigTypeByRef_o}
	${CXX} ${DEUBG_CXX} ${ContigTypeByRef_o} -o ContigTypeByRef

ContigGraphType : clean ${ContigGraphType_o}
	${CXX} ${DEUBG_CXX} ${ContigGraphType_o} -o ContigGraphType

UnicomGraph : clean ${UnicomGraph_o}
	${CXX} ${DEUBG_CXX} ${UnicomGraph_o} -o UnicomGraph

ClusterGap : clean ${ClusterGap_o}
	${CXX} ${DEUBG_CXX} ${ClusterGap_o} -o ClusterGap

job12 : clean ${job_12_o}
	${CXX} ${DEUBG_CXX} ${job_12_o} -o job12

LinearClusterResult : clean ${LinearClusterResult_o}
	${CXX} ${DEUBG_CXX} ${LinearClusterResult_o} -o LinearClusterResult

MergeClusterResult : clean ${MergeClusterResult_o}
	${CXX} ${DEUBG_CXX} ${MergeClusterResult_o} -o MergeClusterResult


LinkCheck: clean ${LinkCheck_o}
	${CXX} ${DEUBG_CXX} ${LinkCheck_o} -o LinkCheck

LinkCheck_1: clean ${LinkCheck_1_o}
	${CXX} ${DEUBG_CXX} ${LinkCheck_1_o} -o LinkCheck_1


ContigRoadFill2Fasta : clean  ${ContigRoadFill2Fasta_o}
	${CXX} ${DEUBG_CXX} ${ContigRoadFill2Fasta_o} -o ContigRoadFill2Fasta

SolveMulti : clean ${SolveMulti_o}
	${CXX} ${DEUBG_CXX} ${SolveMulti_o} -o SolveMulti

ClusterLengthSim : clean ${ClusterLengthSim_o}
	${CXX} ${DEUBG_CXX} ${ClusterLengthSim_o} -o ClusterLengthSim

dirty	   =${jobs_o} ${jobs}

clean:
	rm -rf ${dirty}
