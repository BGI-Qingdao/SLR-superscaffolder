#include "biocommon/fastqa/fasta.h"


namespace BGIQD{
namespace FASTQA{

    SeqItemFactory::Process SeqItemFactory::process( const std::string &line, BGIQD::FASTQA::SeqItem &item )
    {
        if( type == Type::FASTA )
        {
            return process_fasta(line,item);
        }
        else if ( type == Type::FASTQ )
        {
            return process_fastq(line , item );
        }
        return Process::UNKNOW;
    }

    SeqItemFactory::Process SeqItemFactory::process_fasta( const std::string &line, BGIQD::FASTQA::SeqItem &item )
    {
        if( line.empty() )
            return curr;
        Process next1 ;
        if( line[0] == '>' )
        {
            next1 = Line1 ;
        }
        if( curr == Process::UNKNOW )
        {
            curr = Line1 ;
            item.name = 
        }
    }
    SeqItemFactory::Process SeqItemFactory::process_fastq( const std::string &line, BGIQD::FASTQA::SeqItem &item )
    {

    }

} //namespace FASTQA
} //namespace BGIQD 
