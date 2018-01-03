#ifndef __COMMON_SAM_BAM_SAM_PARSER_H__
#define __COMMON_SAM_BAM_SAM_PARSER_H__

#include <string>
#include <vector>

namespace BGIQD{
namespace SAM{

class Head
{
    //TODO 
    int todo;
};

enum CIGAR
{
    NONE = -1,
    M = 0 ,
    I = 1 ,
    D = 2 ,
    N = 3 ,
    S = 4 ,
    H = 5 ,
    P = 6 ,
    EQUAL = 7 ,
    X = 8 ,
};
struct MatchInfo
{
    CIGAR type;
    int start_position_on_read;
    int end_position_on_read;
    int start_position_on_ref;
    int end_position_on_ref;
};

struct MatchDetail
{
    std::vector<MatchInfo> infos;
};

struct MatchData
{
    std::string read_name;
    int flag ;
    std::string ref_name;
    size_t first_match_position;
    int  quality;
    MatchDetail detail;
    int read_len;
    //TODO : other columns . 
};// class MapData
/******************************************************************************
 *
 * Parse 1 raw sam data into struct.
 *
 * ***************************************************************************/
class LineParser
{
    public:
        LineParser( const std::string & line ) : m_line(line) {}
        bool IsVaid() const { return m_line.size() > 0 ; }
        bool IsHead() const { return m_line.at(0) == '@' ; }
        MatchData ParseAsMatchData() const ;
        //TODO 
        Head ParseAsHead() const { return Head() ; }
    private:
        const std::string m_line;
    private:
        size_t ParseStringAsCIGAR(const std::string & str,size_t first_match_on_ref, MatchDetail & detail) const ;
};// class LineParse


} // namespace SAM
} // namespace BGIQD
#endif // __COMMON_SAM_BAM_SAM_PARSER_H__
