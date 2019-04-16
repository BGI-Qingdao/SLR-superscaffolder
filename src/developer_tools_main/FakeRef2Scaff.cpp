#include "biocommon/fasta/fasta.h"
#include <sstream>
#include <string>

bool checkArgs( int argc , char **argv )
{
    if( argc > 1 )
    {
        std::string argv1(argv[1]);
        if( argv1  == "h" 
                ||  argv1  == "-h"
                ||  argv1  == "help"
                ||  argv1  == "--help" )
        {
            std::cerr<<"Usage : \n\t"
                <<argv[0]
                <<" <ref.fa >xxx.sim.scaff"
                <<std::endl;
            return false ;
        }
        else
        {
            std::cerr<<"ERROR : argument is not needed !! "<<std::endl;
            std::cerr<<"Usage : \n\t"
                <<argv[0]
                <<" <ref.fa >xxx.sim.scaff"
                <<std::endl;
            return false ;
        }
    }
    return true ;
}

typedef BGIQD::FASTA::Id_Desc_Head  Header ;
typedef BGIQD::FASTA::Fasta<Header> Fa;

struct RefContig
{
    int start_pos ; /* 1 base */
    int end_pos ;   /* 1 base */

    void Init() 
    {
        start_pos = 0 ;
        end_pos = 0 ;
    }
    bool Valid() const 
    {
        return end_pos > start_pos && start_pos > 0 ;
    }
    bool Empty() const 
    {
        return start_pos == 0 && end_pos == 0 ;
    }
};

struct SimScaff
{
    private:
        int lstart ;  /* 1 base */
        int lend ;    /* 1 base */
        int rstart ;  /* 1 base */
        int rend ;    /* 1 base */

        void Init() 
        {
            lstart = 0 ;
            lend = 0 ;
            rstart = 0 ;
            rend = 0;
        }

        int nlen() const { return rstart-lend-1; }

    public:

        bool Valid()  const
        {
            return rend > rstart 
                && rstart > lend+1  
                && lend >lstart
                && lstart > 0 ;
        }

        Fa ToReal( const Fa & ref , int index ) const 
        {
            assert(nlen() > 0 );
            Fa ret ;
            std::ostringstream ost;
            ost <<">scaffold "<<index
                <<"\tL:"<<lstart<<','<<lend
                <<"\tN:"<<nlen()
                <<"\tR:"<<rstart<<','<<rend
                ;
            ret.AddHead( ost.str() );
            ret.AddSeq(ref.seq.atcgs.substr(lstart-1,lend-lstart+1));
            ret.AddSeq(std::string( nlen() , 'N' ) );
            ret.AddSeq(ref.seq.atcgs.substr(rstart-1 , rend-rstart + 1 ));
            return ret ;
        }

        static SimScaff GenSimScaff( 
                const RefContig & ref_contig ,
                int start_pos ,
                int & next_start_pos 
                )
        {
            const int left_len = 2000 ;
            const int right_len = 2000 ;
            const int n_len_max = 1000 ;

            int need_space = left_len + right_len + n_len_max ;
            SimScaff ret ;
            ret.Init() ;
            if( ref_contig.end_pos - start_pos +1 < need_space  )
                return ret ;
            ret.lstart = start_pos ;
            ret.lend = start_pos + left_len -1 ;
            ret.rstart = ret.lend + 1 + (std::rand() % n_len_max +1);
            ret.rend = ret.rstart + right_len -1 ;
            next_start_pos = ret.rend +1 ;
            return ret ;
        }
};

int main(int argc , char **argv)
{
    // Check parameters
    if( !checkArgs( argc , argv ) )
    {
        return -1 ;
    }

    // Load ref 
    std::vector<Fa> AllRef;
    BGIQD::FASTA::FastaReader<Fa> Reader ;
    Reader.LoadAllFasta(std::cin ,AllRef); 

    if( AllRef.size() < 1 )
    {
        std::cerr<<"no ref sequence loaded !!! error !!! exit ... "<<std::endl;
        return -1 ;
    }

    if( AllRef.size() > 1 )
    {
        std::cerr<<"more than 1 ref sequence loaded !!! error !!! exit ... "<<std::endl;
        std::cerr<<"this tool only suppert 1 ref sequence!!!"<<std::endl;
        return -1 ;
    }
    // Split ref by N
    const auto & the_ref =  AllRef[0] ;
    const auto &seq = the_ref.seq.atcgs ;
    bool n_prev = false ;
    int curr_pos = 0 ;
    std::vector<RefContig> AllRefContig;
    RefContig tmp ;
    tmp.Init() ;
    for( int i = 0 ; i < (int)seq.size() ; i ++ )
    {
        curr_pos ++ ;
        char x = seq.at(i);
        if( x == 'n' || x == 'N' )
        {
            if( ! n_prev )
            {
                if( tmp.Valid() )
                {
                    AllRefContig.push_back(tmp);
                    tmp.Init();
                }
                n_prev = true ;
            }
        }
        else
        {
            if( n_prev )
            {
                assert( tmp.Empty() );
                tmp.start_pos = curr_pos ;
                n_prev = false ;
            }
            else
                tmp.end_pos = curr_pos ;
        }
    }
    if( tmp.Valid() )
        AllRefContig.push_back(tmp);

    // Gen AllScaff
    std::vector<SimScaff> scaffs;
    for( const auto & a_ref_contig : AllRefContig )
    {
        int start_pos = a_ref_contig.start_pos ;
        while(1)
        {
            auto a_scaff = SimScaff::GenSimScaff(a_ref_contig,start_pos,start_pos);
            if( a_scaff.Valid() )
                scaffs.push_back(a_scaff);
            else
                break ;
        }
    }
    // Print scaff
    int index = 0 ;
    for( const auto & a_scaff : scaffs )
    {
        index ++ ;
        auto scaff_fa = a_scaff.ToReal(the_ref,index);
        std::cout<<scaff_fa.head.Head()<<std::endl;
        std::cout<<scaff_fa.seq.Seq(100);
    }

    return 0;
}
