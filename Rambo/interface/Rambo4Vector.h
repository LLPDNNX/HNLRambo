#ifndef HNLRambo_Rambo_Rambo4Vector_h
#define HNLRambo_Rambo_Rambo4Vector_h

#include "TLorentzVector.h"

#include <array>

class Rambo4Vector
{
    protected:
        std::array<double,4> epxpypz_;
    
    public:
        Rambo4Vector():
            epxpypz_{{0,0,0,0}}
        {
        }
        
        inline double& operator[](size_t i)
        {
            return epxpypz_[i];
        }
        
        inline double operator [](size_t i) const
        {
            return epxpypz_[i];
        }
        
        TLorentzVector toLorentzVector() const
        {
            TLorentzVector vector;
            vector.SetPxPyPzE(
                epxpypz_[1],
                epxpypz_[2],
                epxpypz_[3],
                epxpypz_[0]
            );
            return vector;
        }
};

#endif
