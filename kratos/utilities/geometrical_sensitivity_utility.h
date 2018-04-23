//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED)
#define KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED

// System includes
#include <memory>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) GeometricalSensitivityUtility
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(GeometricalSensitivityUtility);

    typedef boost::numeric::ublas::matrix<double> MatrixType;

    typedef MatrixType JacobianType;

    typedef MatrixType ShapeFunctionsLocalGradientType;

    typedef MatrixType ShapeFunctionsGradientType;
    
    typedef unsigned IndexType;

    typedef boost::numeric::ublas::indirect_array<DenseVector<std::size_t>> IndirectArrayType;

    typedef boost::numeric::ublas::matrix_indirect<const MatrixType, IndirectArrayType> SubMatrixType;

    template <class T>
    using matrix_row = boost::numeric::ublas::matrix_row<T>;

    ///@}
    ///@name Life Cycle
    ///@{

    GeometricalSensitivityUtility(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De);

    ///@}
    ///@name Operations
    ///@{

    void CalculateSensitivity(IndexType iNode, IndexType iCoord, double& rDetJ_Deriv, ShapeFunctionsGradientType& rDN_DX_Deriv) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    const JacobianType& mrJ;
    const ShapeFunctionsLocalGradientType& mrDN_De;
    MatrixType mCofactorJ;
    double mDetJ;

    ///@}
    ///@name Private Operations
    ///@{

    void Initialize();

    double CalculateDeterminantOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    MatrixType CalculateCofactorOfJacobianSensitivity(IndexType iNode, IndexType iCoord) const;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* KRATOS_GEOMETRICAL_SENSITIVITY_UTILITY_H_INCLUDED defined */
