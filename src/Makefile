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
				-I./stLFR\
				-lz\

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
				common/args/argsparser.cpp\

test_cpp   =	test.cpp\
				common/time/timetools_test.cpp\
				common/string/stringtools_test.cpp\
				biocommon/sam_bam/sam_parser_test.cpp\
				common/args/argsparser_test.cpp\

job_01_cpp =	job01.cpp
job_01_o   =	${job_01_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_03_cpp =	job03.cpp
job_03_o   =	${job_03_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_02_cpp =	job02.cpp
job_02_o   =	${job_02_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_04_cpp =	job04.cpp
job_04_o   =	${job_04_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_05_cpp =	job05.cpp
job_05_o   =	${job_05_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_06_cpp =	job06.cpp
job_06_o   =	${job_06_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_07_cpp =	job07.cpp
job_07_o   =	${job_07_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_08_cpp =	job08.cpp
job_08_o   =	${job_08_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_09_cpp =	job09.cpp
job_09_o   =	${job_09_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_10_cpp =	job10.cpp
job_10_o   =	${job_10_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_11_cpp =	job11.cpp
job_11_o   =	${job_11_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_12_cpp =	job12.cpp
job_12_o   =	${job_12_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_13_cpp =	job13.cpp
job_13_o   =	${job_13_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_14_cpp =	job14.cpp
job_14_o   =	${job_14_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


Sam2ReadInContig_cpp 	= Sam2ReadInContig.cpp
Sam2ReadInContig_o 	 	= ${Sam2ReadInContig_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

test_o 	   =	${test_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

jobs 	   =	job01 job02 job03 job04 job05 job06 job07 job08 job09 job10 \
				job11 job12 job13 job14 Sam2ReadInContig
jobs_o 	   =	${job_01_o} ${job_03_o} ${job_02_o}  ${job_04_o} ${job_05_o}\
				${job_06_o}  ${job_07_o} ${job_08_o} ${job_09_o}  ${job_10_o} \
				${job_11_o} ${job_12_o} ${job_13_o} ${job_14_o} ${Sam2ReadInContig_o}

all : test ${jobs}


job01 : clean ${job_01_o}
	${CXX} ${DEUBG_CXX} ${job_01_o} -o job01
	
job02 : clean ${job_02_o}
	${CXX} ${DEUBG_CXX} ${job_02_o} -o job02 

job03 : clean ${job_03_o}
	${CXX} ${DEUBG_CXX} ${job_03_o} -o job03

job04 : clean ${job_04_o}
	${CXX} ${DEUBG_CXX} ${job_04_o} -o job04

job05 : clean ${job_05_o}
	${CXX} ${DEUBG_CXX} ${job_05_o} -o job05

job06 : clean ${job_06_o}
	${CXX} ${DEUBG_CXX} ${job_06_o} -o job06

job07 : clean ${job_07_o}
	${CXX} ${DEUBG_CXX} ${job_07_o} -o job07

job08 : clean ${job_08_o}
	${CXX} ${DEUBG_CXX} ${job_08_o} -o job08

job09 : clean ${job_09_o}
	${CXX} ${DEUBG_CXX} ${job_09_o} -o job09

job10 : clean ${job_10_o}
	${CXX} ${DEUBG_CXX} ${job_10_o} -o job10

job11 : clean ${job_11_o}
	${CXX} ${DEUBG_CXX} ${job_11_o} -o job11

job12 : clean ${job_12_o}
	${CXX} ${DEUBG_CXX} ${job_12_o} -o job12

job13 : clean ${job_13_o}
	${CXX} ${DEUBG_CXX} ${job_13_o} -o job13

job14 : clean ${job_14_o}
	${CXX} ${DEUBG_CXX} ${job_14_o} -o job14

Sam2ReadInContig: clean ${Sam2ReadInContig_o}
	${CXX} ${DEUBG_CXX} ${Sam2ReadInContig_o} -o Sam2ReadInContig 

test: clean ${test_o} ${source_h}
	${CXX} ${DEUBG_CXX} ${test_o} -o test

dirty	   =	${test_o} ${jobs_o} ${jobs}  test

clean:
	rm -rf ${dirty}
