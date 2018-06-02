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

#if !defined(NORMAL_CDF_H_INCLUDED )
#define  NORMAL_CDF_H_INCLUDED


// System includes
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup MeshingApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Normal Cumulative Distribution Function
/**
 * Implementation of a double precision approximation
 * for the normal CDF. See:
 * Graeme West, 2005, Better approximations to cumulative normal functions
 */
class NormalCDF
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NormalCDF
    KRATOS_CLASS_POINTER_DEFINITION(NormalCDF);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NormalCDF(double Mean, double StDev);

    /// Destructor.
    virtual ~NormalCDF();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    double Probability(double x);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static constexpr double mSqrt2Pi = 2.5066282746310002; //std::sqrt(4.0*std::acos(0.0));

    static constexpr double mSplit = 7.07106781186547;

    static constexpr double mN0 = 220.206867912376;
    static constexpr double mN1 = 221.213596169931;
    static constexpr double mN2 = 112.079291497871;
    static constexpr double mN3 = 33.912866078383;
    static constexpr double mN4 = 6.37396220353165;
    static constexpr double mN5 = 0.700383064443688;
    static constexpr double mN6 = 3.52624965998911e-02;

    static constexpr double mD0 = 440.413735824752;
    static constexpr double mD1 = 793.826512519948;
    static constexpr double mD2 = 637.333633378831;
    static constexpr double mD3 = 296.564248779674;
    static constexpr double mD4 = 86.7807322029461;
    static constexpr double mD5 = 16.064177579207;
    static constexpr double mD6 = 1.75566716318264;
    static constexpr double mD7 = 8.83883476483184e-02;

    ///@}
    ///@name Member Variables
    ///@{

    double mMean;
    double mStDev;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NormalCDF& operator=(NormalCDF const& rOther);

    /// Copy constructor.
    NormalCDF(NormalCDF const& rOther);

    ///@}

}; // Class NormalCDF

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                NormalCDF& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const NormalCDF& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // NORMAL_CDF_H_INCLUDED  defined
