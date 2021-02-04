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
