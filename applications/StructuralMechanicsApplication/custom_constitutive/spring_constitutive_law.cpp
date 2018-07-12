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
#include <pybind11/pybind11.h>
#include <pybind11/eval.h>

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/spring_constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/**
 * Flags related to the CL computation
 */
// Avoiding using the macro since this has a template parameter. If there was no template plase use the KRATOS_CREATE_LOCAL_FLAG macro
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_MASS(Kratos::Flags::Create(0));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_INERTIA_X(Kratos::Flags::Create(1));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_INERTIA_Y(Kratos::Flags::Create(2));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_INERTIA_Z(Kratos::Flags::Create(3));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_STIFFNESS_X(Kratos::Flags::Create(4));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_STIFFNESS_Y(Kratos::Flags::Create(5));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_STIFFNESS_Z(Kratos::Flags::Create(6));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATTIONAL_STIFFNESS_X(Kratos::Flags::Create(7));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATTIONAL_STIFFNESS_Y(Kratos::Flags::Create(8));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATTIONAL_STIFFNESS_Z(Kratos::Flags::Create(9));

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
        "nodal_mass"                 : null,
        "nodal_inertia"              : [null, null, null],
        "nodal_stiffness"            : [null, null, null],
        "nodal_rotational_stiffness" : [null, null, null],
        "interval"                   : [0.0, 1e30]
    })" );

    NewParameters.ValidateAndAssignDefaults(default_parameters);

    /// Setting the flags
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_MASS, NewParameters["nodal_mass"].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_X, NewParameters["nodal_inertia"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_Y, NewParameters["nodal_inertia"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_Z, NewParameters["nodal_inertia"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_X, NewParameters["nodal_stiffness"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_Y, NewParameters["nodal_stiffness"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_Z, NewParameters["nodal_stiffness"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_X, NewParameters["nodal_rotational_stiffness"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_Y, NewParameters["nodal_rotational_stiffness"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_Z, NewParameters["nodal_rotational_stiffness"][2].IsNull());

    /// Creating the functions // TODO: Define this, use pointers so save memory

    /// Getting intervals
    mTimeInterval[0] = NewParameters["interval"][0].GetDouble();
    mTimeInterval[1] = NewParameters["interval"][1].GetDouble();
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
        return mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_MASS);
    }

    return false;
}
/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    // TODO: Define in case we want more complex behaviours
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<Matrix>& rThisVariable)
{
    // TODO: Define in case we want more complex behaviours
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool SpringConstitutiveLaw<TDim>::Has(const Variable<array_1d<double, 3>>& rThisVariable)
{
    if (rThisVariable == NODAL_STIFFNESS) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if  (rThisVariable == NODAL_INERTIA) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if (rThisVariable == NODAL_ROTATIONAL_STIFFNESS) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATTIONAL_STIFFNESS_Z) ) {
            return false;
        } else {
            return true;
        }
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
