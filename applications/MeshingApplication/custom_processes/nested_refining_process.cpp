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


// External includes


// Project includes

#include "custom_processes/nested_refining_process.h"

namespace Kratos {

NestedRefiningProcess::NestedRefiningProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrModelPart(rThisModelPart),
    mParameters(ThisParameters)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "error_estimation_variables"   : [],
        "step_refining_level"          : 2,
        "refining_boundary_condition"  : "Condition2D3N",
        "echo_level"                   : 3
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mEchoLevel = mParameters["echo_level"].GetInt();
}


}  // namespace Kratos.
