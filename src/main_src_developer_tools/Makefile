.PHONY: all clean dev_bin

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
				../biocommon/paf/PAF.cpp\
                ../biocommon/align_common/align_result.cpp\
                ../biocommon/agp/agp.cpp\
		   		../common/files/file_writer.cpp \
		   		../common/files/gzstream.cpp \
				../common/log/log.cpp\
				../common/log/logfilter.cpp\
				../common/time/timetools.cpp\
				../common/string/stringtools.cpp\
				../stLFR/readName2Barcode.cpp\
				../stLFR/barcodeOnContig.cpp\
				../stLFR/ContigCluster.cpp\
                ../stLFR/ScaffInfo.cpp\
				../stLFR/LineGroup.cpp\
				../stLFR/ContigOverlap.cpp\
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

jobs = barcode2scaffold BarcodeDepth bin_ref cpSpeed FakeRef2Scaff gapfillerIDY GenSentieonScaff GenBinSimMatrix IDYtool IDYtool_paf Cluster2MST2DOT Points4GapSize

all :  ${jobs}




barcode2scaffold_cpp = barcode2scaffold.cpp
barcode2scaffold_o = barcode2scaffold.o
barcode2scaffold : clean ${barcode2scaffold_o} ${source_o} ../dev_bin
	${CXX} ${barcode2scaffold_o} ${source_o} ${DEUBG_CXX} -o barcode2scaffold
	mv $@ ../dev_bin/




BarcodeDepth_cpp = BarcodeDepth.cpp
BarcodeDepth_o = BarcodeDepth.o
BarcodeDepth : clean ${BarcodeDepth_o} ${source_o} ../dev_bin
	${CXX} ${BarcodeDepth_o} ${source_o} ${DEUBG_CXX} -o BarcodeDepth
	mv $@ ../dev_bin/




bin_ref_cpp = bin_ref.cpp
bin_ref_o = bin_ref.o
bin_ref : clean ${bin_ref_o} ${source_o} ../dev_bin
	${CXX} ${bin_ref_o} ${source_o} ${DEUBG_CXX} -o bin_ref
	mv $@ ../dev_bin/




cpSpeed_cpp = cpSpeed.cpp
cpSpeed_o = cpSpeed.o
cpSpeed : clean ${cpSpeed_o} ${source_o} ../dev_bin
	${CXX} ${cpSpeed_o} ${source_o} ${DEUBG_CXX} -o cpSpeed
	mv $@ ../dev_bin/




FakeRef2Scaff_cpp = FakeRef2Scaff.cpp
FakeRef2Scaff_o = FakeRef2Scaff.o
FakeRef2Scaff : clean ${FakeRef2Scaff_o} ${source_o} ../dev_bin
	${CXX} ${FakeRef2Scaff_o} ${source_o} ${DEUBG_CXX} -o FakeRef2Scaff
	mv $@ ../dev_bin/




gapfillerIDY_cpp = gapfillerIDY.cpp
gapfillerIDY_o = gapfillerIDY.o
gapfillerIDY : clean ${gapfillerIDY_o} ${source_o} ../dev_bin
	${CXX} ${gapfillerIDY_o} ${source_o} ${DEUBG_CXX} -o gapfillerIDY
	mv $@ ../dev_bin/




GenSentieonScaff_cpp = GenSentieonScaff.cpp
GenSentieonScaff_o = GenSentieonScaff.o
GenSentieonScaff : clean ${GenSentieonScaff_o} ${source_o} ../dev_bin
	${CXX} ${GenSentieonScaff_o} ${source_o} ${DEUBG_CXX} -o GenSentieonScaff
	mv $@ ../dev_bin/




GenBinSimMatrix_cpp = GenBinSimMatrix.cpp
GenBinSimMatrix_o = GenBinSimMatrix.o
GenBinSimMatrix : clean ${GenBinSimMatrix_o} ${source_o} ../dev_bin
	${CXX} ${GenBinSimMatrix_o} ${source_o} ${DEUBG_CXX} -o GenBinSimMatrix
	mv $@ ../dev_bin/




IDYtool_cpp = IDYtool.cpp
IDYtool_o = IDYtool.o
IDYtool : clean ${IDYtool_o} ${source_o} ../dev_bin
	${CXX} ${IDYtool_o} ${source_o} ${DEUBG_CXX} -o IDYtool
	mv $@ ../dev_bin/




IDYtool_paf_cpp = IDYtool_paf.cpp
IDYtool_paf_o = IDYtool_paf.o
IDYtool_paf : clean ${IDYtool_paf_o} ${source_o} ../dev_bin
	${CXX} ${IDYtool_paf_o} ${source_o} ${DEUBG_CXX} -o IDYtool_paf
	mv $@ ../dev_bin/




Cluster2MST2DOT_cpp = Cluster2MST2DOT.cpp
Cluster2MST2DOT_o = Cluster2MST2DOT.o
Cluster2MST2DOT : clean ${Cluster2MST2DOT_o} ${source_o} ../dev_bin
	${CXX} ${Cluster2MST2DOT_o} ${source_o} ${DEUBG_CXX} -o Cluster2MST2DOT
	mv $@ ../dev_bin/




Points4GapSize_cpp = Points4GapSize.cpp
Points4GapSize_o = Points4GapSize.o
Points4GapSize : clean ${Points4GapSize_o} ${source_o} ../dev_bin
	${CXX} ${Points4GapSize_o} ${source_o} ${DEUBG_CXX} -o Points4GapSize
	mv $@ ../dev_bin/


jobs_o=  ${barcode2scaffold_o} ${BarcodeDepth_o} ${bin_ref_o} ${cpSpeed_o} ${FakeRef2Scaff_o} ${gapfillerIDY_o} ${GenSentieonScaff_o} ${GenBinSimMatrix_o} ${IDYtool_o} ${IDYtool_paf_o} ${Cluster2MST2DOT_o} ${Points4GapSize_o}

dirty	   =${jobs_o} ${jobs} ${source_o}

../dev_bin:
	mkdir -p ../dev_bin

clean:
	rm -rf ${dirty}

