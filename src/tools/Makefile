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
				../soap2/contigGraph.cpp\
				../soap2/contigFasta.cpp\
				../soap2/contigType.cpp\

source_o		= ${source_cpp:%.cpp=%.o}

stLFR_2_10X_cpp 	=	stLFR_2_10X.cpp
stLFR_2_10X_o   =	${stLFR_2_10X_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

cpSpeed_cpp 	=	cpSpeed.cpp
cpSpeed_o   =	${cpSpeed_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

trunk_mark_cpp 	=	trunk_mark.cpp
trunk_mark_o   =	${trunk_mark_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

seedGapFreq_cpp 	=	seedGapFreq.cpp
seedGapFreq_o   =	${seedGapFreq_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

gap_cluster_check_cpp 	=	gap_cluster_check.cpp
gap_cluster_check_o   =	${gap_cluster_check_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	stLFR_2_10X cpSpeed trunk_mark\
				gap_cluster_check\
				seedGapFreq

jobs_o 	   =	${stLFR_2_10X_o} ${cpSpeed_o}\
				${gap_cluster_check_o}\
				${seedGapFreq_o}

all : ${jobs}
	
gap_cluster_check : clean ${gap_cluster_check_o}
	${CXX} ${DEUBG_CXX} ${gap_cluster_check_o} -o gap_cluster_check

seedGapFreq : clean ${seedGapFreq_o}
	${CXX} ${DEUBG_CXX} ${seedGapFreq_o} -o seedGapFreq

trunk_mark : clean ${trunk_mark_o}
	${CXX} ${DEUBG_CXX} ${trunk_mark_o} -o trunk_mark

cpSpeed : clean ${cpSpeed_o}
	${CXX} ${DEUBG_CXX} ${cpSpeed_o} -o cpSpeed

stLFR_2_10X : clean ${stLFR_2_10X_o}
	${CXX} ${DEUBG_CXX} ${stLFR_2_10X_o} -o stLFR_2_10X

dirty	   =${jobs_o} ${jobs} 

clean:
	rm -rf ${dirty}
