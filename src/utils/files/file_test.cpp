#include "utils/unittest/Test.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"

TEST_MODULE_INIT(FILES)

TEST(ReaderUnExisitFile)
{
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName("nowayexist.xyzacvbasdasdad");
    CHECK(true , (NULL ==in ) );
    if( in )
        delete in;
}

TEST(WriteThenRead)
{
    auto out= BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName("temp.test.txt");
    CHECK(true, (NULL != out) );
    (*out)<<"Hello world"<<std::endl;
    delete out; out=NULL;

    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName("temp.test.txt");
    CHECK(true , (NULL !=in ) );
    auto check_hello = [](const std::string & line) {
        CHECK("Hello world",line);
    };
    BGIQD::FILES::FileReaderFactory::EachLine(*in,check_hello);
    delete in; in = NULL;
}

TEST(GZIP_WriteThenRead)
{
    auto out= BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName("temp.test.gz");
    CHECK(true, (NULL != out) );
    (*out)<<"Hello world"<<std::endl;
    delete out; out=NULL;

    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName("temp.test.gz");
    CHECK(true , (NULL !=in ) );
    auto check_hello = [](const std::string & line) {
        CHECK("Hello world",line);
    };
    BGIQD::FILES::FileReaderFactory::EachLine(*in,check_hello);
    delete in; in = NULL;
}
