#ifndef HNLRambo_Rambo_rambo_h
#define HNLRambo_Rambo_rambo_h

#include <vector>
#include <random>
#include <array>
#include <exception>

using namespace std;

class Rambo
{
    public:
        typedef std::mt19937 RandomEngine;
    protected:
        RandomEngine rndEngine_;
        std::uniform_real_distribution<double> dist_;
        
        inline double rnd()
        {
            return dist_(rndEngine_);
        }
        
    public:
        Rambo(size_t seed=12345);
        std::vector<std::array<double,4>> generate(
            double et, 
            const std::vector<double>& xm,
            double& wt
        );
};
#endif
