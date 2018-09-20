#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

int main()
{
    int curr_size = 0 ;
    int count = 1 ;
    std::string line , cache ;
    std::cerr<<"Edge_num n n"<<'\n';
    std::cerr<<"index   length  reverseComplement"<<'\n';

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
            for(int i = 0 ; i <(int) cache.size() ;)
            {
                std::cout<<cache[i];
                i++;
                if( i %100 == 0 && i != (int)cache.size() )
                {
                    std::cout<<'\n';
                }
            }
            curr_size = 0 ;
        }

        curr_size = 0 ;
        cache = "";
    }
    return 0;
}

