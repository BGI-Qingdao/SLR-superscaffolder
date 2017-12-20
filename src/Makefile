.PHONY: all clean debug test release

CC 		   =	gcc
CXX 	   =	g++

CXXFLAGS   =	-I./\
				-I./common/sam_bam\
				-I./common/files\
				-I./common/log\
				-I./common/time\
				-I./common/test\
				-lz\
				-g

source_h   =	contig_barcode.h  \
		   		common/files/file_reader.h \
		   		common/files/file_writer.h  \
		   		common/files/gzstream.h \
		   		common/sam_bam/sam_parser.h\
				common/log/log.h\
				common/log/logfilter.h\
				common/time/timetools.h\
				common/test/Test.h\
				common/test/Check.h\

source_cpp =	contig_barcode.cpp \
		   		common/files/file_reader.cpp \
		   		common/files/file_writer.cpp \
		   		common/files/gzstream.cpp \
		   		common/sam_bam/sam_parser.cpp\
				common/log/log.cpp\
				common/log/logfilter.cpp\
				common/time/timetools.cpp

test_cpp   =	test.cpp\
				common/time/timetools_test.cpp\
				common/sam_bam/sam_parser_test.cpp

main_cpp   =	main.cpp

main_o 	   =	${main_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}

test_o 	   =	${test_cpp:%.cpp=%.o} ${source_cpp:%.cpp=%.o}


.cpp.o:
	${CXX} ${CXXFLAGS} -c $< -o $@


all: ${main_o} ${source_h}

test: clean ${test_o} ${source_h}
	${CXX} ${CXXFLAGS} ${test_o} -o test

dirty	   =	${test_o} ${main_o} 

clean:
	rm -rf ${dirty}
