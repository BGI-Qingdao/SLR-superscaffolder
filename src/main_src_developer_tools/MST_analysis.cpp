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

struct MST_AnalysisEdge {
    unsigned int from ;
    unsigned int to ;
    int rank ; // 0 means non-unique relative edge
};

struct MST_AnalysisNode {

    unsigned int id ;

    bool is_unique ;

    enum GraphType {
        UnknowNode = 0 ,
        LinearNode = 1 ,
        TipNode = 2 ,
        JunctionNode = 3 ,
    };
    GraphType graph_type ;
    enum JunctionType {
        UnknowJunction = 0 ,
        TipJunction = 1,
        BrankJunction = 2 ,
        MixedJunction = 3 ,
        NeibJunction = 4 
    };

    JunctionType junction_type ;

    std::vector<int> edges ; // only store the index of edges .

    int BranchNum () const { return edges.size() ; }

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
std::map<unsigned int , MST_AnalysisNode> nodes;
std::vector<MST_AnalysisEdge> edges;
std::map<unsigned int , std::map<unsigned int , float > > JS_matrix;
//       ref                       pos   contig
std::map<std::string , std::vector<std::tuple<int , unsigned int> > > ref_ranks ;
/********************************************************
  *
  * logic functions
  *
  ******************************************************/

void parse_mintree( std::istream & ist )
{

}

/*
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
    if( index + step >= (int)ref_ranks[ref].size() ) return 0 ;
    return std::get<1>(ref_ranks[ref][index+step]);
}

float JS( unsigned int left ,unsigned int r )
{
    if( JS_matrix.find(left) != JS_matrix.end() && JS_matrix[left].find(r) != JS_matrix[left].end() )
        return JS_matrix[left][r];
    else
        return -1 ;
}
// count and print rank "step" 's simularity
int processRank(int step)
{
    unsigned int  ret = 0 ;
    std::string file = output_dir+"/step_"+std::to_string(step)+"_sim.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & pair : ref_ranks ) 
    {
        for( int i =step ; i< (int)pair.second.size() ; i++ )
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
// count and print rank > "step" 's simularity
int ProcessRankMoreThan(int step)
{
    unsigned int  ret = 0 ;
    std::string file = output_dir+"/step_morethan_"+std::to_string(step)+"_sim.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        for( const auto & pair1 : pair.second ) 
        {
            unsigned int right = pair1.first;
            if( left > right ) 
                continue ;
            if( contigs.find( left ) == contigs.end() )
                continue ;
            if( contigs.find( right ) == contigs.end() )
                continue ;
            const auto & left_contig = contigs.at(left);
            const auto & right_contig = contigs.at(right);
            if( left_contig.ref == "" || right_contig.ref =="" ||  left_contig.ref != right_contig.ref )
                continue ;
            if( abs( left_contig.rank - right_contig.rank ) <= step )
                continue ;
            if ( left_contig.rank > right_contig.rank )
                std::swap(left,right);
            {
                ret ++ ;
                (*out)<<left<<'\t'<<right<<'\t'<<pair1.second<<'\n';
            }
        }
    }
    delete out;
    return ret ;
}
// count and print cross ref simularity
int ProcessCrossRef()
{
    unsigned int  ret = 0 ;
    std::string file = output_dir+"/cross_ref_sim.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        for( const auto & pair1 : pair.second ) 
        {
            unsigned int right = pair1.first ;
            if( left > right ) 
                continue ;
            if( contigs.find( left ) == contigs.end() )
                continue ;
            if( contigs.find( right ) == contigs.end() )
                continue ;
            const auto & left_contig = contigs.at(left);
            const auto & right_contig = contigs.at(right);
            if( left_contig.ref == "" || right_contig.ref =="" ||  left_contig.ref == right_contig.ref )
                continue ;
            {
                ret ++ ;
                (*out)<<left<<'\t'<<right<<'\t'<<pair1.second<<'\n';
            }
        }
    }
    delete out;
    return ret ;
}

// count and print cross ref simularity
int ProcessNonSeeds()
{
    unsigned int  ret = 0 ;
    std::string file = output_dir+"/non_seeds_sim.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        for( const auto & pair1 : pair.second ) 
        {
            unsigned int right = pair1.first ;
            if( left > right ) 
                continue ;
            if( contigs.find( left ) == contigs.end() )
                continue ;
            if( contigs.find( right ) == contigs.end() )
                continue ;
            const auto & left_contig = contigs.at(left);
            const auto & right_contig = contigs.at(right);
            if( !(left_contig.ref == "" ||  right_contig.ref =="") )
                continue ;
            {
                ret ++ ;
                (*out)<<left<<'\t'<<right<<'\t'<<pair1.second<<'\n';
            }
        }
    }
    delete out;
    return ret ;
}

void PrintSortSims(int step)
{
    std::map< int , std::vector<float> > step_sim_map;
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        for( const auto & pair1 : pair.second ) 
        {
            unsigned int right = pair1.first ;
            if( left > right ) 
                continue ;
            if( contigs.find( left ) == contigs.end() )
                continue ;
            if( contigs.find( right ) == contigs.end() )
                continue ;
            const auto & left_contig = contigs.at(left);
            const auto & right_contig = contigs.at(right);
            if( left_contig.ref == "" ||  right_contig.ref =="" )
            {
                step_sim_map[-1].push_back(pair1.second);
            }
            else if ( left_contig.ref != right_contig.ref )
            {
                step_sim_map[-2].push_back(pair1.second);
            }
            else {
                int curr_step = std::abs( left_contig.rank - right_contig.rank );
                if( curr_step > step )
                {
                    step_sim_map[-3].push_back(pair1.second);
                }
                else
                {
                    step_sim_map[curr_step].push_back(pair1.second);
                }
            }
        }
    }
    for( auto & pair : step_sim_map )
    {
        std::sort( pair.second.rbegin() , pair.second.rend() );
    }
    std::string file = output_dir+"/sorted_sim.csv";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    int total_column = 0 ;;
    // print header
     for( const auto & pair : step_sim_map ) {
        if( pair.first == -3 ) {
            (*out)<<"rank > "<<step<<" , ";
        }
        else if ( pair.first == -2 ) {
            (*out)<<"corss_ref ,";
        }
        else if ( pair.first == -1 ) {
            (*out)<<"non_unique ,";
        }
        else {
            (*out)<<"rank "<<pair.first<<" ,";
        }
        total_column ++ ;
    }
    (*out)<<'\n';
    // print data matrix
    int end_column = 0 ;
    size_t index = 0 ;
    while( end_column < total_column )
    {
        end_column = 0 ;
        for( const auto & pair : step_sim_map )
        {
            if( pair.second.size() > index )
            {
                (*out)<<pair.second[index];
            }
            else end_column ++ ;
            (*out)<<" , ";
        }
        index ++ ;
        (*out)<<"\n";
    }
    delete out ;
}

int printMaxSimInfo(int step)
{
    int ret = 0 ;
    auto get_max_sims = []( const std::map<unsigned int , float> & sims ) {
        assert(sims.size() > 0 ) ;
        unsigned int right = -1 ; float max_sim = 0 ;
        for( const auto & pair : sims )
        {
            if( pair.second > max_sim ){
                right = pair.first ;
                max_sim = pair.second ;
            }
        }
        return std::make_pair(right,max_sim);
    };
    std::string file = output_dir+"/max_sim_by_rank_"+std::to_string(step)+"_.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    (*out)<<"contig_id\tJS\tparter_contig_id\n";
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        auto pair1 = get_max_sims(pair.second);
        unsigned int right = pair1.first ;
        if( contigs.find( left ) == contigs.end() )
            continue ;
        if( contigs.find( right ) == contigs.end() )
            continue ;
        const auto & left_contig = contigs.at(left);
        const auto & right_contig = contigs.at(right);
        if( left_contig.ref == "" ||  right_contig.ref =="" )
            continue ;
        else if ( left_contig.ref != right_contig.ref )
            continue ;
        else {
            int curr_step = std::abs( left_contig.rank - right_contig.rank );
            if( curr_step != step )
                continue ;
            ret ++ ;
            (*out)<<left<<'\t'<<pair1.second<<'\t'<<right<<'\n';
        }
    }
    delete out ;
    return ret ;
}

int printMaxSimInfoNot12()
{
    int ret = 0 ;
    auto get_max_sims = []( const std::map<unsigned int , float> & sims ) {
        assert(sims.size() > 0 ) ;
        unsigned int right = -1 ; float max_sim = 0 ;
        for( const auto & pair : sims )
        {
            if( pair.second > max_sim ){
                right = pair.first ;
                max_sim = pair.second ;
            }
        }
        return std::make_pair(right,max_sim);
    };
    std::string file = output_dir+"/max_sim_by_rank_not_1_2_.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    (*out)<<"contig_id\tJS\tparter_contig_id\n";
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        auto pair1 = get_max_sims(pair.second);
        unsigned int right = pair1.first ;
        if( contigs.find( left ) == contigs.end() )
            continue ;
        if( contigs.find( right ) == contigs.end() )
            continue ;
        const auto & left_contig = contigs.at(left);
        const auto & right_contig = contigs.at(right);
        if( left_contig.ref == "" ||  right_contig.ref =="" )
        {
            ret ++ ;
            (*out)<<left<<'\t'<<pair1.second<<'\t'<<right<<'\n';
            continue ;
        }
        else if ( left_contig.ref != right_contig.ref )
        {
            ret ++ ;
            (*out)<<left<<'\t'<<pair1.second<<'\t'<<right<<'\n';
            continue ;
        }
        else {
            int curr_step = std::abs( left_contig.rank - right_contig.rank );
            if( curr_step <= 2 )
                continue ;
            ret ++ ;
            (*out)<<left<<'\t'<<pair1.second<<'\t'<<right<<'\n';
        }
    }
    delete out ;
    return ret ;
}

void PrintContigSims( int step )
{
    std::map< int , std::map< int , std::vector<float> > > contig_step_sim_map;
    for( const auto & pair : JS_matrix )
    {
        unsigned int left = pair.first;
        for( const auto & pair1 : pair.second ) 
        {
            unsigned int right = pair1.first ;
            if( left > right )
                continue ;
            if( contigs.find( left ) == contigs.end() )
                continue ;
            if( contigs.find( right ) == contigs.end() )
                continue ;
            const auto & left_contig = contigs.at(left);
            const auto & right_contig = contigs.at(right);
            if( left_contig.ref == "" ||  right_contig.ref =="" )
            {
                contig_step_sim_map[left][-1] .push_back( pair1.second );
                contig_step_sim_map[right][-1] .push_back( pair1.second );
            }
            else if ( left_contig.ref != right_contig.ref )
            {
                contig_step_sim_map[left][-2] .push_back(  pair1.second );
                contig_step_sim_map[right][-2] .push_back( pair1.second );
            }
            else {
                int curr_step = std::abs( left_contig.rank - right_contig.rank );
                if( curr_step > step )
                {
                    contig_step_sim_map[left][-3] .push_back( pair1.second );
                    contig_step_sim_map[right][-3] .push_back( pair1.second );
                }
                else
                {
                    contig_step_sim_map[left][curr_step] .push_back( pair1.second );
                    contig_step_sim_map[right][curr_step] .push_back( pair1.second );
                }
            }
        }
    }
    //auto get_rank1 = [](const std::map< int , std::vector<float> > & buff , int step  ) -> float {
    //    if( buff.find(step) == buff.end() ) return 0 ;
    //    float ret = 0 ;
    //    for( float x : buff.at(step) )
    //        if( x > ret ) ret = x ;
    //    return ret ;
    //};
    //std::vector< std::tuple< float ,int > > index_buffer;
    //for( const auto & pair : contig_step_sim_map )
    //{
    //    index_buffer.push_back( std::make_tuple( get_rank1(pair.second , 1 ) , pair.first ) ) ;
    //}
    //std::sort( index_buffer.rbegin() , index_buffer.rend() );

    std::string file = output_dir+"/contig_id_sim.csv";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    //int total_column = 3 + step ;
    // print header
    (*out)<<"contig ,";
    (*out)<<"rank > "<<step<<" , ";
    (*out)<<"corss_ref ,";
    (*out)<<"non_unique ,";
    for( int i =0 ; i <step ; i++ )
        (*out)<<"rank "<<i+1<<" ,";
    (*out)<<'\n';
    // print data matrix
    //int end_column = 0 ;
    //size_t index = 0 ;
    for( const auto & pair : contig_step_sim_map )
    {
        int contig = pair.first ;
        for( const auto & pair1 : pair.second  )
        {
            int type = pair1.first ;
            for( float sim : pair1.second ) 
            {
                (*out)<<contig<<" , ";
                if( type == -3 ){
                    (*out)<<sim<<" , , , ";
                    for( int i = 0 ; i < step ; i ++ ) (*out) << ", ";
                } else if ( type == -2 ) {
                    (*out)<<" , "<<sim<<" , , ";
                    for( int i = 0 ; i < step ; i ++ ) (*out) << ", ";
                } else if ( type == -1 ) {
                    (*out)<<" , , "<<sim<<" , ";
                    for( int i = 0 ; i < step ; i ++ ) (*out) << ", ";
                } else {
                    (*out)<<" , , ,  ";
                    for( int i = 0 ; i < step ; i ++ ) {
                        if( i == type -1 ) (*out) << sim <<"," ; else (*out) << ", ";
                    }
                }
                (*out)<<'\n';
            }
        }
    }
    delete out ;
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
*/

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
/*
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
*/


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
        DEFINE_ARG_OPTIONAL(int , step_max , "will print detail about step [1,step_max] ","3");
    END_PARSE_ARGS

    if( test.to_bool() )
    {
        //test_sorted_unique_contigs();
        return 0;
    }
    check_file_read(log_mst.to_string());
    check_file_read(contig_index.to_string());
    check_file_read(mst_seeds.to_string());
    check_file_read(sorted_unique.to_string());
    check_file_read(mst_cluster.to_string());
    create_output_dir(output_dir.to_string() );
    /*
    report("############    Summary report  #####################\n");
    { // log_mst
        auto mst = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(log_mst.to_string());
        report("# contig-sim sub graph  :   "+std::to_string(parse_log_mst(*mst)));
        delete mst ;
    }
    { // xxx.ContigIndex 
        auto index = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(contig_index.to_string());
        report("# contigs               :   "+std::to_string(parse_contig_index(*index)));
        delete index ;
    }
    { // xxx.mst.seeds
        auto seeds = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(mst_seeds.to_string());
        report("# seed-contigs          :   "+std::to_string(parse_mst_seeds(*seeds)));
        delete seeds ;
    }
    { // sort_unique_contig.txt
        auto sorts = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(sorted_unique.to_string());
        report("# unique-contigs        :   "+std::to_string(parse_sort_unique_contigs(*sorts)));
        delete sorts ;
    }
        report("# unique-seeds          :   "+std::to_string(count_seeds_uniq()));
        report("# non-unque-seeds       :   "+std::to_string(count_seeds_wrong()));
    prepare_ref_ranks();
    {//  xxx.mst.cluster
        auto cluster = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(mst_cluster.to_string());
        parse_mst_cluster(*cluster);
        delete cluster ;
    }
    */
    report("############    Summary end    #####################\n");
    report("END");
    return 0 ;
}
