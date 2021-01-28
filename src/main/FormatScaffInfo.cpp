#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/misc/Error.h"

#include "stLFR/ScaffInfo.h"

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_OPTIONAL(bool, format_index ,"format the scaff_index row & scaff_id column","no");
        DEFINE_ARG_OPTIONAL(bool, format_startpos ,"format the start_pos row & scaff_id column","no");
    END_PARSE_ARGS;
    BGIQD::stLFR::ScaffInfoHelper helper;
    helper.LoadAllScaff(std::cin);

    if ( format_index.to_bool() )
    {
        helper.FormatAllIndex();
    }
    if ( format_startpos.to_bool() )
    {
        helper.FormatAllStartPos();
    }
    helper.PrintAllScaff(std::cout);
    return 0;
}
