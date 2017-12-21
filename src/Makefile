.PHONY: all clean 

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-I./\
				-I./biocommon/sam_bam\
				-I./common/files\
				-I./common/log\
				-I./common/time\
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

source_cpp =	contig_barcode.cpp \
		   		common/files/file_reader.cpp \
		   		common/files/file_writer.cpp \
		   		common/files/gzstream.cpp \
		   		biocommon/sam_bam/sam_parser.cpp\
				common/log/log.cpp\
				common/log/logfilter.cpp\
				common/time/timetools.cpp\
				common/string/stringtools.cpp\

test_cpp   =	test.cpp\
				common/time/timetools_test.cpp\
				common/string/stringtools_test.cpp\
				biocommon/sam_bam/sam_parser_test.cpp

main_cpp   =	main.cpp

main_o 	   =	${main_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

test_o 	   =	${test_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${DEUBG_CXX} -c $< -o $@

all : test debug release

debug : clean ${main_o}
	${CXX} ${DEUBG_CXX} ${main_o} -o debug 
	
release : clean ${main_o}
	${CXX} ${RELEASE_CXX} ${main_o} -o release

test: clean ${test_o} ${source_h}
	${CXX} ${DEUBG_CXX} ${test_o} -o test

dirty	   =	${test_o} ${main_o} 

clean:
	rm -rf ${dirty}
