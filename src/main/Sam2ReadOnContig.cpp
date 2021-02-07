/**********************************************************
 *
 * @Brief  :
 *
 *   Parse SAM and print in customized format.
 *
 * *******************************************************/
#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"
#include "utils/multithread/MultiThread.h"
#include "utils/misc/contigIndex.h"
#include "utils/misc/fileName.h"
#include "utils/misc/TagId.h"
#include "stLFR/EasySam.h"

#include "stLFR/stLFRRead.h"
#include <iostream>
#include <string>
#include <cassert>

// SAM header. ignored most of the time.
struct Head
{
    enum HeadType
    {
        Unknow = -1,
        HeadLine= 0 ,
        Sequence = 1 ,
        ReadGroup = 2 ,
        Program = 3,
        OneLineComment = 4
    };

    struct SequenceData
    {
        std::string name ;
        int length ;
    };

    struct VersionData
    {
        std::string version;
    };
    struct Data
    {
        SequenceData sequenceData;
        VersionData  versionData;
    };
    HeadType type ;
    Data d ;
};

// The CIGAR string
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

// One piece of cigar like 10M or 3I or 4D.
struct MatchInfo
{
    CIGAR type;
    int start_position_on_read;  // this index start from 0
    int end_position_on_read;
    int start_position_on_ref;   // this index start from 1
    int end_position_on_ref;
    int len ;
    // start and end construct a [ start , end ] area .
};

// flags of SAM
struct FLAGS
{
    int ox1:1;
    int ox2:1;
    int ox4:1;
    int ox8:1;
    int ox10:1;
    int ox20:1;
    int ox40:1;
    int ox80:1;
    int ox100:1;
    int ox200:1;
    int ox400:1;
    int ox800:1;
};

// take the advantage of union to parse flag
union FLAGS_INT
{
    FLAGS flags;
    int num;
};

// A CIGAR will be parsed as one MatchDetail
struct MatchDetail
{
    std::vector<MatchInfo> infos;
};

// One alignment!
struct MatchData
{
    std::string read_name;
    FLAGS_INT flags;
    //int flag ;
    std::string ref_name;
    size_t first_match_position;
    int  quality;
    MatchDetail detail;
    int read_len;
    bool origin ;
    std::string next_ref_name ;
    int next_ref_pos;
    int insert_size ;
    bool XA;
    bool MD;
    MatchData() : read_name("") , ref_name("")
                  , first_match_position(0)
                  , quality(0) , read_len(-1){}

    int CalcRead1Position() const;
    bool IsP() const ;
    bool IsE() const ;
    bool IsPrimaryMatch() const ;
    bool IsPCRduplicae() const ;
    bool IsPEInSameRef() const ;
    bool IsPEBothMatch() const ;
    bool IsPEBothProperlyMatch() const ;
    bool IsSupplementaryMatch() const ;
    bool IsSecondaryMatch() const ;
    bool Valid() const { return ! detail.infos.empty() ; }
    bool UnMap() const ;
    bool OtherUnMap() const ;
    bool IsReverseComplete() const ;
    int firstMatchInRefNoReverse() const ;

    // the match len + insert len + delete len 
    int total_result_len() const ;

    // the match len + insert len + delete len 
    int total_match_len() const ;

    int total_in_len()  const ;

    int total_del_len() const ;

    int total_clip_len() const ;

    int total_indel_len() const ;
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
        Head ParseAsHead() const ;

    private:
        const std::string m_line;
    private:
        size_t ParseStringAsCIGAR(const std::string & str,size_t first_match_on_ref, MatchDetail & detail) const ;
};// class LineParse



int MatchData::total_del_len()  const 
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::D 
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}

int MatchData::total_in_len()  const 
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::I 
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}


int MatchData::total_clip_len()  const 
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::S
                || x.type == CIGAR::H 
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}
int MatchData::total_indel_len()  const 
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::I 
                || x.type == CIGAR::D 
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}

int MatchData::total_result_len() const
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::M 
                || x.type == CIGAR::X
                || x.type == CIGAR::EQUAL
                || x.type == CIGAR::D 
                || x.type == CIGAR::I 
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}

int MatchData::total_match_len() const
{
    int ret = 0 ;
    for ( const  auto & x : detail.infos )
    {
        if ( x.type == CIGAR::M 
                || x.type == CIGAR::X
                || x.type == CIGAR::EQUAL
           )
        {
            ret += x.len;
        }
    }
    return ret ;
}

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

bool MatchData::IsSecondaryMatch() const 
{
    return ((flags.flags.ox100 & 0x1) == 1 && (flags.flags.ox800 & 0x1) == 0 );
}
bool MatchData::IsSupplementaryMatch() const 
{
    return ((flags.flags.ox800 & 0x1) == 1 );
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
            if( detail.infos[i].type == CIGAR::M )
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
                return (int)info.start_position_on_ref + read_len - (int)info.start_position_on_read - 1;
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
    data.XA = false ;
    data.MD = false ;
    while( ! ist.eof() )
    {
        ist>>extra;
        if(extra.size() >1 &&  extra[0]  == 'M' && extra[1] == 'D' ) 
        {
            data.MD = true ;
        }
        if(extra.size() >1 &&  extra[0]  == 'X' && extra[1] == 'A' )
        {
            data.XA = true ;
        }

        extra.clear();
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
            info_buffer.len = 0 ;
            if ( read_move > 0 )
            {
                info_buffer.start_position_on_read = curr_position_on_read;
                info_buffer.end_position_on_read = curr_position_on_read + read_move -1 ;
                ret = info_buffer.end_position_on_read + 1;
                curr_position_on_read += read_move ;
                info_buffer.len =  read_move ;
            }
            if ( ref_move > 0 )
            {
                info_buffer.start_position_on_ref = curr_position_on_ref ;
                info_buffer.end_position_on_ref = curr_position_on_ref + ref_move -1 ;
                curr_position_on_ref += ref_move ;
                info_buffer.len =  ref_move;
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
                    append_new_info( CIGAR::D, -1, length );
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


//
// Struct to wrap all global variables and functions
//
struct AppConfig
{
    BGIQD::LOG::logger loger;
    std::string barcode_2_num_file;
    BGIQD::MISC::FileNames fName;
    bool has_barcode_in_read_name ;

    typedef BGIQD::stLFR::stLFRHeader stLFRHeader ;

    void Init(const std::string & prefix  ,const std::string & b2n_f, bool b )
    {
        // init loger
        loger.Init("Sam2ReadInContig");
        barcode_2_num_file = b2n_f ;
        fName.Init(prefix);
        has_barcode_in_read_name = ! b;
    }

    BGIQD::MISC::StringIdCache barcodeIds;

    BGIQD::MISC::StringIdCache readNameIds ;

    void LoadBarcode2Num()
    {
        BGIQD::LOG::timer t (loger,"LoadBarcode2Num");
        barcodeIds.preload = true ;
        barcodeIds.Load(fName.barcodeList());
    }

    void LoadRead2Num()
    {
        BGIQD::LOG::timer t(loger,"LoadRead2Num");
        readNameIds.preload = true ;
        readNameIds.Load(fName.readNameList());
    }

    void ParseSam2ReadOnContig()
    {
        std::vector<BGIQD::EASY_SAM::EasySam> easy_cache;
        auto sam_in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig_sam());
        // basic function
        // long long readId = 1 ;
        auto print1read= [&](const MatchData &d)
        {
            BGIQD::EASY_SAM::EasySam tmp;
            //tmp.read_id = readId++;
            tmp.contig_name = std::stoul(d.ref_name);
            tmp.pos_1bp = d.CalcRead1Position();
            tmp.match_reverse = d.IsReverseComplete() ;
            stLFRHeader header;
            header.Init(d.read_name);
            assert(header.type != stLFRHeader::ReadType::Unknow );
            tmp.barcode = barcodeIds.Id(header.barcode_str);
            tmp.read_id = readNameIds.Id(header.readName);
            tmp.is_p = d.IsP();
            if( tmp.is_p == 0 )
                tmp.read_id ++ ;
            tmp.pe_match = (d.IsPEBothMatch() && ! d.XA);
            tmp.insert_size = d.insert_size ;
            easy_cache.push_back(tmp);
        };
        long long count = 0 ;
        auto parseline = [&](const std::string & line)
        {
            LineParser l(line);
            if( ! l.IsVaid() || l.IsHead() )
            {
                return ;
            }
            auto mdata = l.ParseAsMatchData();
            if( mdata.UnMap()) 
                return ;
            if( ! mdata.IsPrimaryMatch() )
                return ;
            count ++ ;
            print1read(mdata);
            if( count % 1000000 == 0 )
                loger<<BGIQD::LOG::lstart()<<count<<"   pair maped reads processed ..."<<BGIQD::LOG::lend();
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*sam_in , parseline);
        delete sam_in ;

        auto b2r_out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.read2contig());
        for( const auto & item : easy_cache)
        {
            (*b2r_out)<<item.ToString()<<'\n';
        }
        delete b2r_out;
    }
}config;

int main(int argc , char ** argv)
{
    // parse args
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,prefix, "prefix. Input xxx.read2contig.sam \n\
                xxx.readNameList\n\
                xxx.barcodeList ;\n\
                Output xxx.read2contig");
    END_PARSE_ARGS

        config.Init(prefix.to_string() , ""/*barcodeList.to_string()*/ ,false );// no_stLFR.to_bool());
    BGIQD::LOG::timer t(config.loger,"Same2ReadOnContig");
    config.LoadRead2Num();
    config.LoadBarcode2Num() ;
    config.ParseSam2ReadOnContig();
    return 0;
}
