#include "biocommon/paf/PAF.h"

#include <sstream>
namespace BGIQD {
    namespace PAF {
        // 15      4432    0       259     +       9125    13577   9       273     242     266     0       tp:A:S  mm:i:15 gn:i:9  go:i:2  cg:Z:136M7D13M2I108M
        void PAF_Item::InitFromString(const std::string &line)
        {
            std::istringstream ist(line);
            ist>>query_name>>query_len>>query_start>>query_end>>query_char;
            ist>>target_name>>target_len>>target_start>>target_end;
            ist>>len_query_match>>len_target_match>>quality;
        }
    }
}
