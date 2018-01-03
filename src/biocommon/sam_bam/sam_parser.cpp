#include "sam_parser.h"
#include <sstream>
namespace BGIQD{
namespace SAM{

MatchData LineParser::ParseAsMatchData() const 
{
    MatchData data;
    std::istringstream ist(m_line);
    std::string cigar;
    ist>>data.read_name>>data.flag>>data.ref_name>>data.first_match_position>>data.quality>>cigar;

    data.read_len = ParseStringAsCIGAR(cigar,data.first_match_position,data.detail);
    return data;
}//ParseAsMatchData

size_t LineParser::ParseStringAsCIGAR( const std::string &str ,size_t first_match_on_ref, MatchDetail & detail) const 
{
    std::string  number_buffer;
    MatchInfo info_buffer;
    size_t curr_position_on_ref = first_match_on_ref;
    size_t curr_position_on_read = 0 ;
    size_t len = str.size() ;
    size_t ret = 0;
    auto append_new_info = [ &info_buffer 
        , &curr_position_on_ref , &curr_position_on_read
        , &detail  , & ret ]
        ( CIGAR type , int read_move  , int ref_move)
        {
            info_buffer.type = type ;
            info_buffer.start_position_on_read = -1 ;
            info_buffer.end_position_on_read = -1 ;
            info_buffer.start_position_on_ref = -1 ;
            info_buffer.end_position_on_ref = -1 ;
            if ( read_move > 0 )
            {
                info_buffer.start_position_on_read = curr_position_on_read;
                info_buffer.end_position_on_read = curr_position_on_read + read_move -1 ;
                ret = info_buffer.end_position_on_read + 1;
                curr_position_on_read += read_move ;
            }
            if ( ref_move > 0 )
            {
                info_buffer.start_position_on_ref = curr_position_on_ref ;
                info_buffer.end_position_on_ref = curr_position_on_ref + ref_move -1 ;
                curr_position_on_ref += ref_move ;
            }

            detail.infos.push_back(info_buffer);
        };

    for(size_t i = 0 ; i < len ; i++ )
    {
        switch(str[i])
        {
            case '0'...'9':
                number_buffer+=str[i];
                break;
            case '*':
                detail.infos.clear();
                append_new_info( CIGAR::NONE , -1 , -1 );
                return ret;
            case 'M':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::M, length , length );
                    number_buffer.clear();
                }
                break;
            case '=':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::EQUAL, length , length );
                    number_buffer.clear();
                }
                break;
            case 'X':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::X, length , length );
                    number_buffer.clear();
                }
                break;
            case 'I':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::I, length , -1);
                    number_buffer.clear();
                }
                break;
            case 'H':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::H, length , -1);
                    number_buffer.clear();
                }
                break;
            case 'S':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::S, length , -1);
                    number_buffer.clear();
                }
                break;
            case 'D':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::I, -1, length );
                    number_buffer.clear();
                }
                break;
            case 'P':
                //TODO :
                number_buffer.clear();
                break;
            case 'N':
                {
                    int length = std::stoi(number_buffer);
                    append_new_info( CIGAR::N, -1 , length );
                    number_buffer.clear();
                }
                break;
        }
    }
    return ret;
}

}//namespace SAM
}//namespace BGIQD
