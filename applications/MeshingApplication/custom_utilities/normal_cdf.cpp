//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes
#include <cmath>


// External includes 


// Project includes
#include "custom_utilities/normal_cdf.h"


namespace Kratos
{

NormalCDF::NormalCDF(double Mean, double StDev) : mMean(Mean), mStDev(StDev){};

double NormalCDF::Probability(double x)
{
    const double z = std::abs(x);
    double c = 0.0;

    if (z <= 37.0)
    {
        const double e = std::exp(-z*z/2.0);
        if (z < mSplit)
        {
            const double n = (((((mN6*z + mN5)*z + mN4)*z + mN3)*z + mN2)*z + mN1)*z + mN0;
            const double d = ((((((mD7*z + mD6)*z + mD5)*z + mD4)*z + mD3)*z + mD2)*z + mD1)*z + mD0;
            c = e*n/d;
        }
        else
        {
            const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
            c = e/(mSqrt2Pi*f);
        }
    }
    return x<=0.0 ? c : 1-c;
}
  
}  // namespace Kratos.


