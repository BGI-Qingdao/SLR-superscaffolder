.PHONY: all clean 

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-std=c++11 -I./\
				-I./biocommon/sam_bam\
				-I./common/files\
				-I./common/log\
				-I./common/time\
				-I./common/args\
				-I./common/test\
				-I./common/string\
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
				common/log/log.cpp\
				common/log/logfilter.cpp\
				common/time/timetools.cpp\
				common/string/stringtools.cpp\
				common/args/argsparser.cpp\

test_cpp   =	test.cpp\
				common/time/timetools_test.cpp\
				common/string/stringtools_test.cpp\
				biocommon/sam_bam/sam_parser_test.cpp\
				common/args/argsparser_test.cpp\

job_01_cpp =	job01.cpp

job_02_cpp =	job02.cpp

job_01_o   =	${job_01_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

job_02_o   =	${job_02_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

test_o 	   =	${test_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

all : test job01 job02

job01 : clean ${job_01_o}
	${CXX} ${DEUBG_CXX} ${job_01_o} -o job01
	
job02 : clean ${job_02_o}
	${CXX} ${DEUBG_CXX} ${job_02_o} -o job02 

test: clean ${test_o} ${source_h}
	${CXX} ${DEUBG_CXX} ${test_o} -o test

dirty	   =	${test_o} ${job_01_o} ${job_02_o} job01 job02 test

clean:
	rm -rf ${dirty}
