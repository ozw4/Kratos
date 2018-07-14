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
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_STIFFNESS_X(Kratos::Flags::Create(7));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_STIFFNESS_Y(Kratos::Flags::Create(8));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_STIFFNESS_Z(Kratos::Flags::Create(9));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_DAMPING_RATIO_X(Kratos::Flags::Create(10));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_DAMPING_RATIO_Y(Kratos::Flags::Create(11));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_DAMPING_RATIO_Z(Kratos::Flags::Create(12));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_DAMPING_RATIO_X(Kratos::Flags::Create(13));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_DAMPING_RATIO_Y(Kratos::Flags::Create(14));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::NULL_ROTATIONAL_DAMPING_RATIO_Z(Kratos::Flags::Create(15));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_MASS(Kratos::Flags::Create(16));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_INERTIA_X(Kratos::Flags::Create(17));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_INERTIA_Y(Kratos::Flags::Create(18));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_INERTIA_Z(Kratos::Flags::Create(19));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_STIFFNESS_X(Kratos::Flags::Create(20));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_STIFFNESS_Y(Kratos::Flags::Create(21));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_STIFFNESS_Z(Kratos::Flags::Create(22));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_STIFFNESS_X(Kratos::Flags::Create(23));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_STIFFNESS_Y(Kratos::Flags::Create(24));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_STIFFNESS_Z(Kratos::Flags::Create(25));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_DAMPING_RATIO_X(Kratos::Flags::Create(26));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_DAMPING_RATIO_Y(Kratos::Flags::Create(27));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_DAMPING_RATIO_Z(Kratos::Flags::Create(28));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_DAMPING_RATIO_X(Kratos::Flags::Create(29));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_DAMPING_RATIO_Y(Kratos::Flags::Create(30));
template<std::size_t TDim>
const Kratos::Flags SpringConstitutiveLaw<TDim>::CONSTANT_ROTATIONAL_DAMPING_RATIO_Z(Kratos::Flags::Create(31));

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
        "nodal_mass"                     : null,
        "nodal_inertia"                  : [null, null, null],
        "nodal_stiffness"                : [null, null, null],
        "nodal_rotational_stiffness"     : [null, null, null],
        "nodal_damping_ratio"            : [null, null, null],
        "nodal_rotational_damping_ratio" : [null, null, null],
        "additional_dependence_variables": [],
        "interval"                       : [0.0, 1e30]
    })" );

    NewParameters.ValidateAndAssignDefaults(default_parameters);

    /// The list of additional additional dependence variables
    const SizeType n_variables = NewParameters["additional_dependence_variables"].size();
    mAdditionalDependenceVariables.resize(n_variables);
    for (IndexType i_var = 0; i_var < n_variables; ++i_var){
        const std::string& variable_name = NewParameters["additional_dependence_variables"].GetArrayItem(i_var).GetString();
        mAdditionalDependenceVariables[i_var] = variable_name;
    }

    /// Setting the flags
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_MASS, NewParameters["nodal_mass"].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_X, NewParameters["nodal_inertia"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_Y, NewParameters["nodal_inertia"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_INERTIA_Z, NewParameters["nodal_inertia"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_X, NewParameters["nodal_stiffness"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_Y, NewParameters["nodal_stiffness"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_STIFFNESS_Z, NewParameters["nodal_stiffness"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_X, NewParameters["nodal_rotational_stiffness"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Y, NewParameters["nodal_rotational_stiffness"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Z, NewParameters["nodal_rotational_stiffness"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_DAMPING_RATIO_X, NewParameters["nodal_damping_ratio"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Y, NewParameters["nodal_damping_ratio"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Z, NewParameters["nodal_damping_ratio"][2].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_X, NewParameters["nodal_rotational_damping_ratio"][0].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Y, NewParameters["nodal_rotational_damping_ratio"][1].IsNull());
    mConstitutiveLawFlags.Set(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Z, NewParameters["nodal_rotational_damping_ratio"][2].IsNull());

    /// Creating the functions
    std::string function_string;
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_MASS)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_MASS, NewParameters["nodal_mass"].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_MASS)) {
            mFunctions.insert({NODAL_MASS.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_mass"].GetString())});
        } else {
            mConstantValues.insert({NODAL_MASS.Key(), NewParameters["nodal_mass"].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_X)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_INERTIA_X, NewParameters["nodal_inertia"][0].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_X)) {
            mFunctions.insert({NODAL_INERTIA_X.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_inertia"][0].GetString())});
        } else {
            mConstantValues.insert({NODAL_INERTIA_X.Key(), NewParameters["nodal_inertia"][0].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_Y)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_INERTIA_Y, NewParameters["nodal_inertia"][1].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_Y)) {
            mFunctions.insert({NODAL_INERTIA_Y.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_inertia"][1].GetString())});
        } else {
            mConstantValues.insert({NODAL_INERTIA_Y.Key(), NewParameters["nodal_inertia"][1].GetDouble()});
        }
    }
    if (TDim == 3) {
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_Z)) {
            mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_INERTIA_Z, NewParameters["nodal_inertia"][2].IsDouble());
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_Z)) {
                mFunctions.insert({NODAL_INERTIA_Z.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_inertia"][2].GetString())});
            } else {
                mConstantValues.insert({NODAL_INERTIA_Z.Key(), NewParameters["nodal_inertia"][2].GetDouble()});
            }
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_X)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_STIFFNESS_X, NewParameters["nodal_stiffness"][0].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_X)) {
            mFunctions.insert({NODAL_STIFFNESS_X.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_stiffness"][0].GetString())});
        } else {
            mConstantValues.insert({NODAL_STIFFNESS_X.Key(), NewParameters["nodal_stiffness"][0].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_Y)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Y, NewParameters["nodal_stiffness"][1].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Y)) {
            mFunctions.insert({NODAL_STIFFNESS_Y.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_stiffness"][1].GetString())});
        } else {
            mConstantValues.insert({NODAL_STIFFNESS_Y.Key(), NewParameters["nodal_stiffness"][1].GetDouble()});
        }
    }
    if (TDim == 3) {
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_Z)) {
            mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Z, NewParameters["nodal_stiffness"][2].IsDouble());
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Z)) {
                mFunctions.insert({NODAL_STIFFNESS_Z.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_stiffness"][2].GetString())});
            } else {
                mConstantValues.insert({NODAL_STIFFNESS_Z.Key(), NewParameters["nodal_stiffness"][2].GetDouble()});
            }
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_X)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_X, NewParameters["nodal_rotational_stiffness"][0].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_X)) {
            mFunctions.insert({NODAL_ROTATIONAL_STIFFNESS_X.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_stiffness"][0].GetString())});
        } else {
            mConstantValues.insert({NODAL_ROTATIONAL_STIFFNESS_X.Key(), NewParameters["nodal_rotational_stiffness"][0].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Y)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Y, NewParameters["nodal_rotational_stiffness"][1].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Y)) {
            mFunctions.insert({NODAL_ROTATIONAL_STIFFNESS_Y.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_stiffness"][1].GetString())});
        } else {
            mConstantValues.insert({NODAL_ROTATIONAL_STIFFNESS_Y.Key(), NewParameters["nodal_rotational_stiffness"][1].GetDouble()});
        }
    }
    if (TDim == 3) {
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Z)) {
            mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Z, NewParameters["nodal_rotational_stiffness"][2].IsDouble());
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Z)) {
                mFunctions.insert({NODAL_ROTATIONAL_STIFFNESS_Z.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_stiffness"][2].GetString())});
            } else {
                mConstantValues.insert({NODAL_ROTATIONAL_STIFFNESS_Z.Key(), NewParameters["nodal_rotational_stiffness"][2].GetDouble()});
            }
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_X)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_X, NewParameters["nodal_damping_ratio"][0].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_X)) {
            mFunctions.insert({NODAL_DAMPING_RATIO_X.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_damping_ratio"][0].GetString())});
        } else {
            mConstantValues.insert({NODAL_DAMPING_RATIO_X.Key(), NewParameters["nodal_damping_ratio"][0].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Y)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Y, NewParameters["nodal_damping_ratio"][1].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Y)) {
            mFunctions.insert({NODAL_DAMPING_RATIO_Y.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_damping_ratio"][1].GetString())});
        } else {
            mConstantValues.insert({NODAL_DAMPING_RATIO_Y.Key(), NewParameters["nodal_damping_ratio"][1].GetDouble()});
        }
    }
    if (TDim == 3) {
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Z)) {
            mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Z, NewParameters["nodal_damping_ratio"][2].IsDouble());
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Z)) {
                mFunctions.insert({NODAL_DAMPING_RATIO_Z.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_damping_ratio"][2].GetString())});
            } else {
                mConstantValues.insert({NODAL_DAMPING_RATIO_Z.Key(), NewParameters["nodal_damping_ratio"][2].GetDouble()});
            }
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_X)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_X, NewParameters["nodal_rotational_damping_ratio"][0].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_X)) {
            mFunctions.insert({NODAL_ROTATIONAL_DAMPING_RATIO_X.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_damping_ratio"][0].GetString())});
        } else {
            mConstantValues.insert({NODAL_ROTATIONAL_DAMPING_RATIO_X.Key(), NewParameters["nodal_rotational_damping_ratio"][0].GetDouble()});
        }
    }
    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Y)) {
        mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Y, NewParameters["nodal_rotational_damping_ratio"][1].IsDouble());
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Y)) {
            mFunctions.insert({NODAL_ROTATIONAL_DAMPING_RATIO_Y.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_damping_ratio"][1].GetString())});
        } else {
            mConstantValues.insert({NODAL_ROTATIONAL_DAMPING_RATIO_Y.Key(), NewParameters["nodal_rotational_damping_ratio"][1].GetDouble()});
        }
    }
    if (TDim == 3) {
        if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Z)) {
            mConstitutiveLawFlags.Set(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Z, NewParameters["nodal_rotational_damping_ratio"][2].IsDouble());
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Z)) {
                mFunctions.insert({NODAL_ROTATIONAL_DAMPING_RATIO_Z.Key(),  Kratos::make_shared<PythonConstitutiveLawFunction>(NewParameters["nodal_rotational_damping_ratio"][2].GetString())});
            } else {
                mConstantValues.insert({NODAL_ROTATIONAL_DAMPING_RATIO_Z.Key(), NewParameters["nodal_rotational_damping_ratio"][2].GetDouble()});
            }
        }
    }

    /// Getting intervals
    mTimeInterval[0] = NewParameters["interval"][0].GetDouble();
    mTimeInterval[1] = NewParameters["interval"][1].GetDouble();
}

//******************************COPY CONSTRUCTOR************************************/
/***********************************************************************************/

template<std::size_t TDim>
SpringConstitutiveLaw<TDim>::SpringConstitutiveLaw(const SpringConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
    , mFunctions(rOther.mFunctions)
    , mConstantValues(rOther.mConstantValues)
    , mNodalIndex(rOther.mNodalIndex)
    , mConstitutiveLawFlags(rOther.mConstitutiveLawFlags)
    , mTimeInterval(rOther.mTimeInterval)
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
bool SpringConstitutiveLaw<TDim>::Has(const Variable<int>& rThisVariable)
{
    if (rThisVariable == NODAL_INDEX) {
        return true;
    }

    return false;
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
    if  (rThisVariable == NODAL_INERTIA) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_INERTIA_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if (rThisVariable == NODAL_STIFFNESS) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_STIFFNESS_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if (rThisVariable == NODAL_ROTATIONAL_STIFFNESS) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if (rThisVariable == NODAL_DAMPING_RATIO) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_DAMPING_RATIO_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Z) ) {
            return false;
        } else {
            return true;
        }
    } else if (rThisVariable == NODAL_ROTATIONAL_DAMPING_RATIO) {
        if (mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_X) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Y) &&
            mConstitutiveLawFlags.Is(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Z) ) {
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
int& SpringConstitutiveLaw<TDim>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    if (rThisVariable == NODAL_INDEX) {
        rValue = mNodalIndex;
    }

    rValue = 0;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void SpringConstitutiveLaw<TDim>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == NODAL_INDEX) {
        mNodalIndex = rValue;
    }
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
    const GeometryType& geom = rParameterValues.GetElementGeometry();
    const ProcessInfoType& process_info = rParameterValues.GetProcessInfo();
    const double time = process_info[TIME];
    if (time >= mTimeInterval[0] && time <= mTimeInterval[1]) {
        if (rThisVariable == NODAL_MASS) {
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_X))
                    rValue = mFunctions[NODAL_MASS.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue = mConstantValues[NODAL_MASS.Key()];
            } else {
                rValue = 0.0;
            }
        } else {
            rValue = 0.0;
        }
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
    if (time >= mTimeInterval[0] && time <= mTimeInterval[1]) {
//     if (Has(rThisVariable)) {
//
//     } else {
        rValue = ZeroVector(2 * TDim);
//     }
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
    if (time >= mTimeInterval[0] && time <= mTimeInterval[1]) {
//     if (Has(rThisVariable)) {
//
//     } else {
        rValue = ZeroMatrix(TDim, TDim);
//     }
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
    const GeometryType& geom = rParameterValues.GetElementGeometry();
    const ProcessInfoType& process_info = rParameterValues.GetProcessInfo();
    const double time = process_info[TIME];
    if (time >= mTimeInterval[0] && time <= mTimeInterval[1]) {
        if (rThisVariable == NODAL_INERTIA) {
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_X))
                    rValue[0] = mFunctions[NODAL_INERTIA_X.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[0] = mConstantValues[NODAL_INERTIA_X.Key()];
            } else
                rValue[0] = 0.0;
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_Y)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_Y))
                    rValue[1] = mFunctions[NODAL_INERTIA_Y.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[1] = mConstantValues[NODAL_INERTIA_Y.Key()];
            } else
                rValue[1] = 0.0;
            if (TDim == 3) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_INERTIA_Z)) {
                    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_INERTIA_Z))
                        rValue[2] = mFunctions[NODAL_INERTIA_Z.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                    else
                        rValue[2] = mConstantValues[NODAL_INERTIA_Z.Key()];
                } else
                    rValue[2] = 0.0;
            } else {
                rValue[2] = 0.0;
            }
        } else if (rThisVariable == NODAL_STIFFNESS) {
           if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_X))
                    rValue[0] = mFunctions[NODAL_STIFFNESS_X.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[0] = mConstantValues[NODAL_STIFFNESS_X.Key()];
           } else
                rValue[0] = 0.0;
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_Y)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Y))
                    rValue[1] = mFunctions[NODAL_STIFFNESS_Y.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[1] = mConstantValues[NODAL_STIFFNESS_Y.Key()];
            } else
                rValue[1] = 0.0;
            if (TDim == 3) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_STIFFNESS_Z)) {
                    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_STIFFNESS_Z))
                        rValue[2] = mFunctions[NODAL_STIFFNESS_Z.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                    else
                        rValue[2] = mConstantValues[NODAL_STIFFNESS_Z.Key()];
                } else
                    rValue[2] = 0.0;
            } else {
                rValue[2] = 0.0;
            }
        } else if (rThisVariable == NODAL_ROTATIONAL_STIFFNESS) {
           if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_X))
                    rValue[0] = mFunctions[NODAL_ROTATIONAL_STIFFNESS_X.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[0] = mConstantValues[NODAL_ROTATIONAL_STIFFNESS_X.Key()];
            } else
                rValue[0] = 0.0;
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Y)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Y))
                    rValue[1] = mFunctions[NODAL_ROTATIONAL_STIFFNESS_Y.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[1] = mConstantValues[NODAL_ROTATIONAL_STIFFNESS_Y.Key()];
            } else
                rValue[1] = 0.0;
            if (TDim == 3) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_STIFFNESS_Z)) {
                    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_STIFFNESS_Z))
                        rValue[2] = mFunctions[NODAL_ROTATIONAL_STIFFNESS_Z.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                    else
                        rValue[2] = mConstantValues[NODAL_ROTATIONAL_STIFFNESS_Z.Key()];
                } else
                    rValue[2] = 0.0;
            } else {
                rValue[2] = 0.0;
            }
        } else if (rThisVariable == NODAL_DAMPING_RATIO) {
           if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_X))
                    rValue[0] = mFunctions[NODAL_DAMPING_RATIO_X.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[0] = mConstantValues[NODAL_DAMPING_RATIO_X.Key()];
            } else
                rValue[0] = 0.0;
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Y)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Y))
                    rValue[1] = mFunctions[NODAL_DAMPING_RATIO_Y.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[1] = mConstantValues[NODAL_DAMPING_RATIO_Y.Key()];
            } else
                rValue[1] = 0.0;
            if (TDim == 3) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_DAMPING_RATIO_Z)) {
                    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_DAMPING_RATIO_Z))
                        rValue[2] = mFunctions[NODAL_DAMPING_RATIO_Z.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                    else
                        rValue[2] = mConstantValues[NODAL_DAMPING_RATIO_Z.Key()];
                } else
                    rValue[2] = 0.0;
            } else {
                rValue[2] = 0.0;
            }
        } else if (rThisVariable == NODAL_ROTATIONAL_DAMPING_RATIO) {
           if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_X)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_X))
                    rValue[0] = mFunctions[NODAL_ROTATIONAL_DAMPING_RATIO_X.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[0] = mConstantValues[NODAL_ROTATIONAL_DAMPING_RATIO_X.Key()];
            } else
                rValue[0] = 0.0;
            if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Y)) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Y))
                    rValue[1] = mFunctions[NODAL_ROTATIONAL_DAMPING_RATIO_Y.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                else
                    rValue[1] = mConstantValues[NODAL_ROTATIONAL_DAMPING_RATIO_Y.Key()];
            } else
                rValue[1] = 0.0;
            if (TDim == 3) {
                if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::NULL_ROTATIONAL_DAMPING_RATIO_Z)) {
                    if (mConstitutiveLawFlags.IsNot(SpringConstitutiveLaw::CONSTANT_ROTATIONAL_DAMPING_RATIO_Z))
                        rValue[2] = mFunctions[NODAL_ROTATIONAL_DAMPING_RATIO_Z.Key()]->CallFunction(geom[mNodalIndex], mAdditionalDependenceVariables, time);
                    else
                        rValue[2] = mConstantValues[NODAL_ROTATIONAL_DAMPING_RATIO_Z.Key()];
                } else
                    rValue[2] = 0.0;
            } else {
                rValue[2] = 0.0;
            }
        } else {
            rValue = ZeroVector(3);
        }
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
    /// Nodal values of interest
    KRATOS_CHECK_VARIABLE_KEY(NODAL_MASS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_INERTIA)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_ROTATIONAL_STIFFNESS)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_DAMPING_RATIO)
    KRATOS_CHECK_VARIABLE_KEY(NODAL_ROTATIONAL_DAMPING_RATIO)

    /// Dependance values
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(ROTATION)
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ANGULAR_ACCELERATION)

    return 0;

}

} // Namespace Kratos
