/**********************************************************
  *
  *
  *
  *
  *
  ********************************************************/


#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "stLFR/CBB.h"

#include <iostream>
#include <sstream>
#include <regex>

#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <map>
#include <vector>
#include <tuple>

/********************************************************
  *
  * basic structures
  *
  ******************************************************/
struct AnalysisContig
{
    unsigned int    id ;
    int             length;
    bool            is_seed ;
    std::string     ref ;
    int             pos ;    // position in ref , the leftest.
    int             rank ;
    bool            is_uniq ;
    float           repeat_factor ;
    bool            ln_1 ;
    bool            ln_2 ;
    bool            ln_3 ;
    void  init() {  length = 0 ; is_seed = false ; pos =0 ; rank = 0; is_uniq = false ; repeat_factor = 0 ; 
                    ln_1 = ln_2 = ln_3 = false ; }
};

struct QuastInfo
{
    int ref_start ;
    int ref_end ;
    int query_start ;
    int query_end ;
    std::string ref ;
    unsigned int contig ;

    void InitFromString(const std::string & line )
    {
        std::string tmp1 ;
        std::istringstream ist(line);
        ist>>ref_start>>ref_end>>query_start>>query_end>>ref>>tmp1;
        std::string ctg ;
        for( auto i : tmp1 )
        {
            if( std::isdigit(i) )
            {
                ctg += i;
            }
            else
                break ;
        }
        contig = std::stoul(ctg);
    }
};

/********************************************************
  *
  * global variables
  *
  ******************************************************/

static std::string output_dir;
std::map<unsigned int , AnalysisContig> contigs ;
std::map<unsigned int , std::map<unsigned int , float > > JS_matrix;
//       ref                       pos   contig
std::map<std::string , std::vector<std::tuple<int , unsigned int> > > ref_ranks ;
/********************************************************
  *
  * logic functions
  *
  ******************************************************/

// @return : sub graph num
int parse_log_mst( std::istream  & ist )
{
    std::string line ;
    std::regex rg(".* into ([0-9]+) sub graph.*");
    int ret ;
    while( ! ist.eof() )
    {
        std::getline(ist,line);
        if(line.empty() ) continue ;
        if( std::regex_match(line , rg) )
        {
            std::smatch rm;
            std::regex_match(line,rm,rg);
            ret=std::stoi(rm[1]);
        }
    }
    return ret ;
}

// @return : total contig number 
int parse_contig_index( std::istream  & ist )
{
    std::string line ;
    std::getline(ist,line);
    std::getline(ist,line);
    int ret = 0 ;
    unsigned int id ; int length;
    while( ! ist.eof() )
    {
        id =0 ; length = 0 ;
        std::getline(ist,line);
        if(line.empty() ) continue ;
        std::istringstream isst(line);
        isst>>id>>length;
        if( length > 0 ) 
        {   assert( contigs.find(id) == contigs.end() );
            contigs[id].init();
            contigs[id].length = length ;
            ret ++ ;
        }
    }
    return ret ;
}
// @return : total seed-contig number 
int parse_mst_seeds( std::istream  & ist )
{
    std::string line ;
    int ret = 0 ;
    unsigned int id ; int length;
    while( ! ist.eof() )
    {
        id =0 ; length = 0 ;
        std::getline(ist,line);
        if(line.empty() ) continue ;
        std::istringstream isst(line);
        isst>>id>>length;
        if( length > 0 && contigs.find(id) != contigs.end() )
        {
            assert( contigs[id].length == length );
            contigs[id].is_seed = true;
            ret ++ ;
        }
    }
    return ret ;
}

int parse_sort_unique_contigs(std::istream & ist)
{
    std::string line ;
    int ret = 0 ;
    while( ! ist.eof() )
    {
        std::getline(ist,line);
        if(line.empty() ) continue ;
        QuastInfo tmp ;
        tmp.InitFromString(line);
        if( contigs.find(tmp.contig) != contigs.end() )
        {
            contigs[tmp.contig].is_uniq = true ;
            contigs[tmp.contig].ref = tmp.ref ;
            contigs[tmp.contig].pos = tmp.ref_start ;
            ret ++ ;
        }
    }
    return ret;
}

int count_seeds_wrong()
{
    int ret = 0 ;
    for( const auto & pair : contigs )
    {
        const auto i = pair.second ;
        if( i.is_seed && ! i.is_uniq ) ret ++ ;
    }
    return ret ;
}
int count_seeds_uniq()
{
    int ret = 0 ;
    for( const auto & pair : contigs )
    {
        const auto i = pair.second ;
        if( i.is_seed && i.is_uniq ) ret ++ ;
    }
    return ret ;
}
void prepare_ref_ranks()
{
    // collect ref data
    for( const auto & pair : contigs )
    {
        const auto i = pair.second ;
        if( i.is_seed && i.is_uniq ) ref_ranks[i.ref].push_back( std::make_tuple( i.pos ,i.id )) ;
    }
    // sort by pos 
    for( auto & pair : ref_ranks )
    {
        std::sort(pair.second.begin() , pair.second.end());
    }
    // re-assign rank to contigs
    for(const auto & pair : ref_ranks )
    {
        int rank = 0 ;
        for( const auto & tup : pair.second )
        {
            contigs[std::get<1>(tup)].rank = rank ;
            rank ++ ;
        }
    }
}

unsigned int nib(const std::string & ref , int index , int step) 
{
    if( index + step < 0 )  return  0 ;
    if( ref_ranks.find( ref ) == ref_ranks.end() ) return 0; 
    if( index + step >= ref_ranks[ref].size() ) return 0 ;
    return std::get<1>(ref_ranks[ref][index+step]);
}

float JS( unsigned int left ,unsigned int r )
{
    if( JS_matrix.find(left) != JS_matrix.end() && JS_matrix[left].find(r) != JS_matrix[left].end() )
        return JS_matrix[left][r];
    else
        return -1 ;
}

int processRank(int step)
{
    unsigned int  ret = 0 ;
    std::string file = output_dir+"/step_"+std::to_string(step)+"_sim.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & pair : ref_ranks ) 
    {
        for( int i =step ; i< pair.second.size() ; i++ )
        {
            int left = std::get<1>(pair.second.at(i-step));
            int r = std::get<1>(pair.second.at(i));
            float sim = JS(left,r);
            if( sim > 0 )
            {
                ret ++ ;
                (*out)<<left<<'\t'<<r<<'\t'<<sim<<'\n';
            }
        }
    }
    delete out;
    return ret ;
}


void parse_mst_cluster(std::istream & ist)
{
    std::string line ;
    while( ! ist.eof() )
    {
        std::getline(ist,line);
        if(line.empty() ) continue ;
        BGIQD::stLFR::ContigRelation tmp;
        tmp.InitFromString(line);
        for( auto x : tmp.sims )
        {
            JS_matrix[tmp.contigId][ x.first      ]= x.second.simularity ;
            JS_matrix[x.first     ][ tmp.contigId ]= x.second.simularity ;
        }
    }
}

void create_output_dir(const std::string & dirname )
{
    if( mkdir ( dirname.c_str() , S_IRUSR|S_IWUSR|S_IXUSR| S_IRGRP | S_IXGRP) != 0  )
    {
        std::cerr<<"FATAL : failed to create directory "<<dirname<<" !\n         exit ... \n";
        exit(1);
    }
    output_dir = dirname ;
}

void report(const std::string &  log)
{
    static std::ostream * report_file = NULL ;
    if( report_file == NULL )
        report_file = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_dir + "/report.txt");
    if( log == "END" && report_file != NULL )
    {
        delete report_file ;
        return ;
    }
    (*report_file) << log <<'\n';
    std::cout<<log<<'\n';
}

/********************************************************
  *
  * check params functions
  *
  ******************************************************/
void check_file_read(const std::string & file_name)
{
    auto t = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file_name);
    if( t == NULL )
    {
        std::cerr<<"FATAL : failed to open "<<file_name<<" to read!"
            <<"\n         exit ... \n";
        exit(1);
    }
    delete t;
}

void check_file_write(const std::string & file_name)
{
    auto t = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file_name);
    if( t == NULL )
    {
        std::cerr<<"FATAL : failed to open "<<file_name<<" to write!"
            <<"\n         exit ... \n";
        exit(1);
    }
    delete t;
}

/********************************************************
  *
  * test function
  *
  ******************************************************/
void test_log_mst()
{
    std::string test="\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       load contig sim graph done\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       contig-sim graph nodes : 2510\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       contig-sim graph edges : 6092\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       split contig sim graph into 43 sub graph\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       CorrectLog      155     155     155     0\n\
       MST     INFO    CST     2109/9/25       16:47:35        :       CorrectLog      257     257     257     0\n";
    std::istringstream ist(test);
    assert( parse_log_mst(ist) == 43 );
    std::cerr<<"test log mst pass!"<<std::endl;
}

void test_contig_index()
{
    std::string test="\n\
        Edge_num n n\n\
        index   length  reverseComplement\n\
        1       5216    1\n\
        3       6305    1\n\
        5       6136    1\n\
        7       3497    1\n\
        9       6023    1\n\
        11      3874    1\n\
        13      2536    1\n\
        15      4432    1\n";
    std::istringstream ist(test);
    assert( parse_contig_index(ist) == 8 );
    std::cerr<<"test contig_index pass!"<<std::endl;
}

void test_mst_seeds()
{
    std::string test="\n\
        1       5216    1\n\
        3       6305    1\n\
        5       6136    1\n\
        7       3497    1\n\
        9       6023    1\n\
        11      3874    1\n\
        15      4432    1\n";
    std::istringstream ist(test);
    assert( parse_mst_seeds(ist) == 7 );
    std::cerr<<"test mst_seeds pass!"<<std::endl;
}

void test_sorted_unique_contigs()
{
    std::string test="\
107988  110691  1       2704    chr19   1_length_2704_cvg_18.0_tip_0        100.0           True\n\
115453  117505  1380    1       chr19   3_length_2556_cvg_18.0_tip_0  96.74           True\n\
242529  246906  1       4387    chr19   7_length_4387_cvg_18.0_tip_0 99.54           True\n\
253077  269825  16752   1       chr19   9_length_16752_cvg_18.0_tip_0       99.51           True\n\
270940  275078  4140    1       chr19   11 99.95           True\n\
274823  282686  7769    1       chr19   15_length_7769_cvg_18.0_tip_0        98.69           True\n\
282487  284843  2356    1       chr19   5_length_2356_cvg_18.0_tip_0 99.96           True\n";
    std::istringstream ist(test);
    assert( parse_sort_unique_contigs(ist) == 7 );
    std::cerr<<"test sort_unique_contigs pass!"<<std::endl;
}

void test_mst_cluster()
{
    std::string test="\
17      87      0.227187        2849    0.367841        7945    0.194999        10989   0.142886\n\
21      891     0.252198        6545    0.133452        8877    0.245732        9035    0.336013        9903    0.403683        103690.203822        11765   0.159417        11785   0.12321\n\
23      2119    0.196319        2189    0.295057        4749    0.44365 5443    0.142857        10611   0.118593\n\
37      111     0.154341        2221    0.157858        2353    0.128563        5031    0.107009        10103   0.339771        102670.257264\n\
43      2621    0.387872        5623    0.111157        5943    0.119642        9067    0.363018\n\
51      5001    0.383734        8755    0.144743        10029   0.220148        10609   0.166975\n\
57      4227    0.440673        7945    0.130128        8909    0.139766        10769   0.233961\n\
65      5121    0.326181        9521    0.404476        9657    0.317382        10231   0.139061        10567   0.425637\n\
69      701     0.120421        1869    0.179641        4733    0.263452        6059    0.387725        11271   0.436948\n";
    std::istringstream ist(test);
    parse_mst_cluster(ist);
    std::cerr<<"test mst_cluster pass!"<<std::endl;
}

/********************************************************
  *
  * main function
  *
  ******************************************************/

int main( int argc , char ** argv )
{

    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,log_mst," log_mst file .");
        DEFINE_ARG_REQUIRED(std::string,contig_index," xxx.ContigIndex file .");
        DEFINE_ARG_REQUIRED(std::string,mst_seeds," xxx.mst.seeds file .");
        DEFINE_ARG_REQUIRED(std::string,sorted_unique," sorted_unique_contigs.txt file .");
        DEFINE_ARG_REQUIRED(std::string,mst_cluster," xxx.mst.cluster file .");
        DEFINE_ARG_REQUIRED(std::string,output_dir," output directoy name .");
        DEFINE_ARG_OPTIONAL(bool, test, "run self test","0");
    END_PARSE_ARGS

    if( test.to_bool() )
    {
        test_log_mst();
        test_contig_index();
        test_mst_seeds();
        test_sorted_unique_contigs();
        return 0;
    }
    check_file_read(log_mst.to_string());
    check_file_read(contig_index.to_string());
    check_file_read(mst_seeds.to_string());
    check_file_read(sorted_unique.to_string());
    check_file_read(mst_cluster.to_string());
    create_output_dir(output_dir.to_string() );

    report("############    Summary report  #####################\n");
    { // log_mst
        auto mst = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(log_mst.to_string());
        report("# contig-sim sub graph  :   "+std::to_string(parse_log_mst(*mst))+"\n");
        delete mst ;
    }
    { // xxx.ContigIndex 
        auto index = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(contig_index.to_string());
        report("# contigs               :   "+std::to_string(parse_contig_index(*index))+"\n");
        delete index ;
    }
    { // xxx.mst.seeds
        auto seeds = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(mst_seeds.to_string());
        report("# seed-contigs          :   "+std::to_string(parse_mst_seeds(*seeds))+"\n");
        delete seeds ;
    }
    { // sort_unique_contig.txt
        auto sorts = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(sorted_unique.to_string());
        report("# unique-contigs        :   "+std::to_string(parse_sort_unique_contigs(*sorts))+"\n");
        delete sorts ;
    }
        report("# unique-seeds          :   "+std::to_string(count_seeds_uniq())+"\n");
        report("# non-unque-seeds       :   "+std::to_string(count_seeds_wrong())+"\n");
    prepare_ref_ranks();
    {//  xxx.mst.cluster
        auto cluster = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(mst_cluster.to_string());
        parse_mst_cluster(*cluster);
        delete cluster ;
    }
        report("# step 1 edge           :   "+std::to_string(processRank(1))+"\n");
        report("# step 2 edge           :   "+std::to_string(processRank(2))+"\n");
        report("# step 3 edge           :   "+std::to_string(processRank(3))+"\n");

    report("############    Summary end    #####################\n");
    report("END");
    return 0 ;
}
