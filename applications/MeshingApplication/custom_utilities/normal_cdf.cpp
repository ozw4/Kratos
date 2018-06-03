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
        if (z < msSplit)
        {
            const double n = (((((msN6*z + msN5)*z + msN4)*z + msN3)*z + msN2)*z + msN1)*z + msN0;
            const double d = ((((((msD7*z + msD6)*z + msD5)*z + msD4)*z + msD3)*z + msD2)*z + msD1)*z + msD0;
            c = e*n/d;
        }
        else
        {
            const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
            c = e/(msSqrt2Pi*f);
        }
    }
    return x<=0.0 ? c : 1-c;
}

}  // namespace Kratos.


