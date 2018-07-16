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
#include "includes/kratos_components.h"
#include "includes/variables.h"
#include "custom_python/python_constitutive_law_function.h"

namespace Kratos
{
namespace Python
{
/**
 * ConstitutiveLawFunction memeber variables definitions
 */
pybind11::object mMainModule;    /// The python main module
pybind11::object mMainNameSpace; /// The variables that generate the dependence of the function
std::string mFunctionBody;       /// The text definting the function

/**
 * ConstitutiveLawFunction methods definitions
 */
ConstitutiveLawFunction::ConstitutiveLawFunction( const std::string& rFunctionBody)
{
    // Compile the function starting from the string function body
    try {
        mMainModule = pybind11::module::import("__main__");
        mMainNameSpace = mMainModule.attr("__dict__");
        pybind11::exec("from math import *", mMainNameSpace);
        mFunctionBody = rFunctionBody;
    } catch(pybind11::error_already_set const&) {
        PyErr_Print();
    }
}

/***********************************************************************************/
/***********************************************************************************/

double ConstitutiveLawFunction::CallFunction(
    const NodeType& ThisNode,
    std::vector<std::string>& AdditionalVariables,
    const double Time
    )
{
    /// Position
    mMainNameSpace["x"] = ThisNode.X();
    mMainNameSpace["y"] = ThisNode.Y();
    mMainNameSpace["z"] = ThisNode.Z();
    mMainNameSpace["X0"] = ThisNode.X0();
    mMainNameSpace["Y0"] = ThisNode.Y0();
    mMainNameSpace["Z0"] = ThisNode.Z0();

    /// Values
    const array_1d<double, 3>& disp = ThisNode.FastGetSolutionStepValue(DISPLACEMENT);
    mMainNameSpace["DISPLACEMENT_X"] = disp[0];
    mMainNameSpace["DISPLACEMENT_Y"] = disp[1];
    mMainNameSpace["DISPLACEMENT_Z"] = disp[2];
    /// Increment of values
    const array_1d<double, 3>& dispDelta = ThisNode.FastGetSolutionStepValue(DISPLACEMENT) - ThisNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
    mMainNameSpace["DISPLACEMENT_DELTA_X"] = dispDelta[0];
    mMainNameSpace["DISPLACEMENT_DELTA_Y"] = dispDelta[1];
    mMainNameSpace["DISPLACEMENT_DELTA_Z"] = dispDelta[2];
    if (ThisNode.SolutionStepsDataHas(VELOCITY_X)) {
        const array_1d<double, 3>& vel = ThisNode.FastGetSolutionStepValue(VELOCITY);
        mMainNameSpace["VELOCITY_X"] = vel[0];
        mMainNameSpace["VELOCITY_Y"] = vel[1];
        mMainNameSpace["VELOCITY_Z"] = vel[2];
        /// Increment of values
        const array_1d<double, 3>& velDelta = ThisNode.FastGetSolutionStepValue(VELOCITY) - ThisNode.FastGetSolutionStepValue(VELOCITY, 1);
        mMainNameSpace["VELOCITY_DELTA_X"] = velDelta[0];
        mMainNameSpace["VELOCITY_DELTA_Y"] = velDelta[1];
        mMainNameSpace["VELOCITY_DELTA_Z"] = velDelta[2];
    }
    if (ThisNode.SolutionStepsDataHas(ACCELERATION_X)) {
        const array_1d<double, 3>& accel = ThisNode.FastGetSolutionStepValue(ACCELERATION);
        mMainNameSpace["ACCELERATION_X"] = accel[0];
        mMainNameSpace["ACCELERATION_Y"] = accel[1];
        mMainNameSpace["ACCELERATION_Z"] = accel[2];
        /// Increment of values
        const array_1d<double, 3>& accelDelta = ThisNode.FastGetSolutionStepValue(ACCELERATION) - ThisNode.FastGetSolutionStepValue(ACCELERATION, 1);
        mMainNameSpace["ACCELERATION_DELTA_X"] = accelDelta[0];
        mMainNameSpace["ACCELERATION_DELTA_Y"] = accelDelta[1];
        mMainNameSpace["ACCELERATION_DELTA_Z"] = accelDelta[2];
    }
    if (ThisNode.SolutionStepsDataHas(ROTATION_X)) {
        const array_1d<double, 3>& theta = ThisNode.FastGetSolutionStepValue(ROTATION);
        mMainNameSpace["ROTATION_X"] = theta[0];
        mMainNameSpace["ROTATION_Y"] = theta[1];
        mMainNameSpace["ROTATION_Z"] = theta[2];
        /// Increment of values
        const array_1d<double, 3>& thetaDelta = ThisNode.FastGetSolutionStepValue(ROTATION) - ThisNode.FastGetSolutionStepValue(ROTATION, 1);
        mMainNameSpace["ROTATION_DELTA_X"] = thetaDelta[0];
        mMainNameSpace["ROTATION_DELTA_Y"] = thetaDelta[1];
        mMainNameSpace["ROTATION_DELTA_Z"] = thetaDelta[2];
    }
    if (ThisNode.SolutionStepsDataHas(ANGULAR_VELOCITY_X)) {
        const array_1d<double, 3>& vang = ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        mMainNameSpace["ANGULAR_VELOCITY_X"] = vang[0];
        mMainNameSpace["ANGULAR_VELOCITY_Y"] = vang[1];
        mMainNameSpace["ANGULAR_VELOCITY_Z"] = vang[2];
        /// Increment of values
        const array_1d<double, 3>& vangDelta = ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY) - ThisNode.FastGetSolutionStepValue(ANGULAR_VELOCITY, 1);
        mMainNameSpace["ANGULAR_VELOCITY_DELTA_X"] = vangDelta[0];
        mMainNameSpace["ANGULAR_VELOCITY_DELTA_Y"] = vangDelta[1];
        mMainNameSpace["ANGULAR_VELOCITY_DELTA_Z"] = vangDelta[2];
    }
    if (ThisNode.SolutionStepsDataHas(ANGULAR_ACCELERATION_X)) {
        const array_1d<double, 3>& aang = ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION);
        mMainNameSpace["ANGULAR_ACCELERATION_X"] = aang[0];
        mMainNameSpace["ANGULAR_ACCELERATION_Y"] = aang[1];
        mMainNameSpace["ANGULAR_ACCELERATION_Z"] = aang[2];
        /// Increment of values
        const array_1d<double, 3>& aangDelta = ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION) - ThisNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 1);
        mMainNameSpace["ANGULAR_ACCELERATION_DELTA_X"] = aangDelta[0];
        mMainNameSpace["ANGULAR_ACCELERATION_DELTA_Y"] = aangDelta[1];
        mMainNameSpace["ANGULAR_ACCELERATION_DELTA_Z"] = aangDelta[2];
    }

    /// Time
    mMainNameSpace["t"] = Time;

    /// Additional variables
    for (auto& variable_name : AdditionalVariables) {
        if (KratosComponents< Variable< array_1d< double, 3> > >::Has(variable_name)) {
            Variable<array_1d< double, 3>> variable = KratosComponents< Variable<array_1d< double, 3>> >::Get(variable_name);
            const array_1d<double, 3>& value_var = ThisNode.FastGetSolutionStepValue(variable);
            std::string aux_string = variable_name + "_X";
            mMainNameSpace[aux_string.c_str()] = value_var[0];
            aux_string = variable_name + "_Y";
            mMainNameSpace[aux_string.c_str()] = value_var[1];
            aux_string = variable_name + "_Z";
            mMainNameSpace[aux_string.c_str()] = value_var[2];
        } else {
            Variable<double> variable = KratosComponents< Variable<double> >::Get(variable_name);
            const double value_var = ThisNode.FastGetSolutionStepValue(variable);
            mMainNameSpace[variable_name.c_str()] = value_var;
        }
    }

    return pybind11::eval(mFunctionBody, mMainNameSpace).cast<double>();
}
} // Namespace Python
} // Namespace Kratos
