#ifndef __BIOCOMMON_FASTQA_FASTA_H__
#define __BIOCOMMON_FASTQA_FASTA_H__
#include <string>
#include <vector>

namespace BGIQD{
namespace FASTQA{

struct SeqItem
{
    std::string name;
    std::string extra_flags;
    std::string seq;
    std::string quality;
};

enum Type
{
    UNKNOW = 0,
    FASTQ = 1,
    FASTA = 2,
} ;

struct SeqDatabase
{
    Type type;
    std::vector<SeqItem> datas;
};

struct SeqItemFactory
{
    enum Process
    {
        UNKNOW = 0 ,
        Line1 = 1 ,
        Line2 = 2 ,
        Line2More = 3 ,
        Line3 = 4,
        Line4 = 5,
        Line4More = 6,
        Finish = 7,
    };

    void Init( Type t) { type = t ; }
    Process process(const std::string & line, SeqItem & item);
    private :
    Type type;
    Process curr;
    SeqItem next;

    Process process_fasta(const std::string & line, SeqItem & item1);
    Process process_fastq(const std::string & line, SeqItem & item1);
};

} //namespace FASTQA
} //namespace BGIQD 

#endif //__BIOCOMMON_FASTQA_FASTA_H__
