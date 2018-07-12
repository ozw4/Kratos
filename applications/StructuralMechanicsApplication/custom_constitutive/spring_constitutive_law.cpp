// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/spring_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

template<std::size_t TDim>
SpringConstitutiveLaw<TDim>::SpringConstitutiveLaw()
    : ConstitutiveLaw()
{
    KRATOS_WARNING("SpringConstitutiveLaw") << "Using default constructor, please use the constructor via parameters" << std::endl;
}

//******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

template<std::size_t TDim>
SpringConstitutiveLaw<TDim>::SpringConstitutiveLaw(Kratos::Parameters NewParameters)
    : ConstitutiveLaw()
{
    Kratos::Parameters default_parameters = Kratos::Parameters(R"(
    {
        "thig"  : "thing"
    })" );

    NewParameters.ValidateAndAssignDefaults(default_parameters);
}

//******************************COPY CONSTRUCTOR************************************/
/***********************************************************************************/

template<std::size_t TDim>
SpringConstitutiveLaw<TDim>::SpringConstitutiveLaw(const SpringConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

//********************************CLONE*********************************************/
/***********************************************************************************/

template<std::size_t TDim>
ConstitutiveLaw::Pointer SpringConstitutiveLaw<TDim>::Clone() const
{
    return Kratos::make_shared<SpringConstitutiveLaw<TDim>>(SpringConstitutiveLaw<TDim>(*this));
}

//*******************************DESTRUCTOR*****************************************/
/***********************************************************************************/

template<std::size_t TDim>
SpringConstitutiveLaw<TDim>::~SpringConstitutiveLaw<TDim>()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/

template<std::size_t TDim>
void SpringConstitutiveLaw<TDim>::GetLawFeatures(Features& rFeatures)
{
    // Set the strain size
    rFeatures.mStrainSize =  2 * TDim;

    // Set the spacedimension
    rFeatures.mSpaceDimension = TDim;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == NODAL_MASS) {
    }

    return false;
}
/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<Vector>& rThisVariable)
{
//     if (rThisVariable == NODAL_MASS) {
//     }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<Matrix>& rThisVariable)
{
//     if (rThisVariable == NODAL_MASS) {
//     }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    if (rThisVariable == NODAL_STIFFNESS) {

    } else if  (rThisVariable == NODAL_INERTIA) {

    } else if (rThisVariable == NODAL_ROTATIONAL_STIFFNESS) {

    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
double& SpringConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (Has(rThisVariable)) {

    } else {
        rValue = 0.0;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Vector& SpringConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (Has(rThisVariable)) {

    } else {
        rValue = ZeroVector(2 * TDim);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Matrix& SpringConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (Has(rThisVariable)) {

    } else {
        rValue = ZeroMatrix(TDim, TDim);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
array_1d<double, 3 > & SpringConstitutiveLaw<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<array_1d<double, 3 > >& rThisVariable,
    array_1d<double, 3 > & rValue
    )
{
    if (Has(rThisVariable)) {

    } else {
        rValue = ZeroVector(3);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
int SpringConstitutiveLaw<TDim>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_INERTIA)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_ROTATIONAL_STIFFNESS)

    return 0;

}

} // Namespace Kratos
