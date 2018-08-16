#include "biocommon/sam_bam/sam_parser.h"
#include <sstream>
#include "common/string/stringtools.h"

namespace BGIQD{
namespace SAM{

bool MatchData::UnMap() const 
{
    return (flags.flags.ox4 & 0x1) == 1 ;
}
bool MatchData::IsPEBothMatch() const 
{
    return next_ref_name != "*";
}
bool MatchData::IsPEInSameRef() const
{
    return next_ref_name == "=" ;
}
bool MatchData::IsPEBothProperlyMatch() const 
{
    return (flags.flags.ox2 & 0x1) == 1 ;
}
bool MatchData::OtherUnMap() const 
{
    return (flags.flags.ox8 & 0x1) == 1 ;
}

bool MatchData::IsP() const
{
    return ( flags.flags.ox40 & 0x1) == 1 ;
}
bool MatchData::IsE() const
{
    return ( flags.flags.ox80 & 0x1 )== 1 ;
}

bool MatchData::IsPrimaryMatch() const
{
    return ((flags.flags.ox800 & 0x1) == 0 ) &&( (flags.flags.ox100 & 0x1 ) == 0 ) ;
}
bool MatchData::IsPCRduplicae() const
{
    return ((flags.flags.ox400 & 0x1) == 1 ) ;
}
bool MatchData::IsReverseComplete() const
{
    return ( flags.flags.ox10 & 0x1 )== 1 ;
}
int MatchData::firstMatchInRefNoReverse() const 
{
    if( ! IsReverseComplete() )
        return (int)first_match_position;
    else
    {
        for(size_t i = 0 ; i < detail.infos.size() ; i++)
        {
            if( detail.infos[i].type == SAM::CIGAR::M )
            {
                return detail.infos[i].start_position_on_ref + read_len - detail.infos[i].start_position_on_read -1  ;
            }
        }
    }
    return -1 ;
}

int MatchData::CalcRead1Position() const
{
    for( int i = 0 ; i < (int) detail.infos.size() ; i ++ )
    {
        const auto & info = detail.infos[i] ;
        if( info.type == CIGAR::EQUAL || info.type == CIGAR::M )
        {
            if( IsReverseComplete() )
            {
                return (int)info.start_position_on_ref + read_len - (int)info.start_position_on_read ;
            }
            else
            {
                return (int)info.start_position_on_ref - (int)info.start_position_on_read ;
            }
        }
    }
    return 0;
}

Head LineParser::ParseAsHead()const
{
    Head h;
    std::istringstream ist(m_line);
    char c;
    ist>>c;//@
    std::string tmp;
    ist>>tmp;
    if ( tmp == "HD" )
    {
        h.type = Head::HeadType::HeadLine;
    }
    else if ( tmp == "SQ")
    {
        h.type = Head::HeadType::Sequence;
        while(1)
        {
            ist>>tmp;
            if(ist.eof())
                break;
            auto t = BGIQD::STRING::split( tmp , ":");

            if( t[0] == "SN" )
            {
                h.d.sequenceData.name = t[1];
            }
            else if (t[0] == "LN")
            {
                h.d.sequenceData.length = std::stoi(t[1]);
            }
        }
    }
    else if ( tmp == "RG" )
    {
        h.type = Head::HeadType::ReadGroup;

    }
    else if ( tmp == "PG" )
    {
        h.type = Head::HeadType::Program;

    }
    else if ( tmp == "CO" )
    {
        h.type = Head::HeadType::OneLineComment;
    }
    return h;
}

MatchData LineParser::ParseAsMatchData() const 
{
    MatchData data;
    std::istringstream ist(m_line);
    std::string cigar;
    ist>>data.read_name
        >>data.flags.num
        >>data.ref_name
        >>data.first_match_position
        >>data.quality
        >>cigar
        >>data.next_ref_name
        >>data.next_ref_pos
        >>data.insert_size
        ;
    std::string extra;
    while( ! ist.eof() )
    {
        ist>>extra;
    }
    if(extra.size() >1 &&  extra[0]  == 'X' && extra[1] == 'A' )
    {
        data.XA = true ;
    }
    else
    {
        data.XA = false ;
    }
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
