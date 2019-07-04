#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
int main()
{
    int curr_size = 0 ;
    int count = 1 ;
    std::string line , cache ;
    std::cerr<<"Edge_num n n"<<'\n';
    std::cerr<<"index   length  reverseComplement"<<'\n';
    std::ofstream ofs;
    ofs.open("fakesoap.name2index.map.txt");
    std::string pre_name ;
    while(! std::getline(std::cin,line).eof() )
    {
        if( line.empty() )
            continue ;
        if( line[0] == '>' )
        {
            if( curr_size > 0)
            {
                std::cout<<">"<<count<<" length "<<curr_size<<" cvg_18.0_tip_0\n";
                std::cerr<<count<<'\t'<<curr_size<<"\t1\n";
                ofs<<pre_name<<'\t'<<count<<'\n';
                count+=2 ;
                for(int i = 0 ; i < (int)cache.size() ;)
                {
                    std::cout<<cache[i];
                    i++;
                    if( i %100 == 0 || i == (int)cache.size() )
                    {
                        std::cout<<'\n';
                    }
                }
                curr_size = 0 ;
                pre_name= "";
            }
            pre_name= "";
            for( int i = 1 ; i <(int)line.size() ; i++ )
            {
                if( std::isblank(line[i] ) )
                        break ;
                pre_name+=line[i] ;
            }
            curr_size = 0 ;
            cache = "";
            continue;
        }
        else
        {
            cache += line;
            curr_size += line.size();
            line="";
        }
    }
    if( curr_size >= 0 )
    {
        if( curr_size > 0)
        {
            std::cout<<">"<<count<<" length "<<curr_size<<" cvg_18.0_tip_0\n";
            std::cerr<<count<<'\t'<<curr_size<<"\t1\n";
            //ofs<<count<<'\t'<<pre_name<<'\n';
            ofs<<pre_name<<'\t'<<count<<'\n';
            for(int i = 0 ; i <(int) cache.size() ;)
            {
                std::cout<<cache[i];
                i++;
                if( i %100 == 0 || i == (int)cache.size() )
                {
                    std::cout<<'\n';
                }
            }
            curr_size = 0 ;
        }

        curr_size = 0 ;
        cache = "";
    }
    ofs.close();
    return 0;
}

