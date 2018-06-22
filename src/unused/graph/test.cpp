#include <iostream>
#include <vector>

int main()
{
    std::vector<int> i;
    i.resize(1000000000);
    std::cout<<i.size()<<std::endl;
    return 0;
}
