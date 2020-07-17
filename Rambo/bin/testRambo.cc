#include "HNLRambo/Rambo/interface/rambo.h"

#include <iostream>

int main()
{
    Rambo rambo(12345);
    
    for (size_t t = 0; t <10; ++t)
    {
        double weight = 0;
        std::vector<std::array<double,4>> momenta = rambo.generate(10.,std::vector<double>{{1.,2.,3.}},weight);
        std::cout<<"w="<<weight<<std::endl;
        
        std::array<double,4> sum{{0,0,0,0}};
        for (size_t i = 0; i < momenta.size(); ++i)
        {
            std::cout<<i<<": ";
            for (size_t k = 0; k < 4; ++k)
            {
                sum[k]+=momenta[i][k];
                std::cout<<momenta[i][k]<<", ";
            }
            std::cout<<" m="<<std::sqrt(
                momenta[i][0]*momenta[i][0]
                -momenta[i][1]*momenta[i][1]
                -momenta[i][2]*momenta[i][2]
                -momenta[i][3]*momenta[i][3]
            )<<std::endl;
        }
        std::cout<<"sum: ";
        for (size_t k = 0; k < 4; ++k)
        {
            std::cout<<sum[k]<<", ";
        }
        std::cout<<" m="<<std::sqrt(
            sum[0]*sum[0]
            -sum[1]*sum[1]
            -sum[2]*sum[2]
            -sum[3]*sum[3]
        )<<std::endl;
    }
    return 0;
};
