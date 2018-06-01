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
#include "processes/fast_transfer_between_model_parts_process.h"
#include "custom_processes/nested_refining_process.h"

namespace Kratos {

NestedRefiningProcess::NestedRefiningProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
) : mrOriginModelPart(rThisModelPart),
    mParameters(ThisParameters)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "own_model_part_name"          : "own",
        "refined_model_part_name"      : "refined",
        "scalar_variable"              : "DISTANCE",
        "gradient_variable"            : "DISTANCE_GRADIENT",
        "probability"                  : 0.9,
        "step_refining_level"          : 2,
        "refining_boundary_condition"  : "Condition2D3N",
        "echo_level"                   : 3
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mEchoLevel = mParameters["echo_level"].GetInt();
    std::string own_model_part_name = mParameters["own_model_part_name"].GetString();
    std::string refined_model_part_name = mParameters["refined_model_part_name"].GetString();

    if (mrOriginModelPart.HasSubModelPart(own_model_part_name))
    {
        mpOwnModelPart = mrOriginModelPart.pGetSubModelPart(own_model_part_name);
    }
    else
    {
        mpOwnModelPart = mrOriginModelPart.CreateSubModelPart(own_model_part_name);
        FastTransferBetweenModelPartsProcess transfer_process = FastTransferBetweenModelPartsProcess(
            *mpOwnModelPart, 
             mrOriginModelPart);
        transfer_process.Execute();
    }

    if (mrOriginModelPart.HasSubModelPart(refined_model_part_name))
    {
        mpRefinedModelPart = mrOriginModelPart.pGetSubModelPart(refined_model_part_name);
    }
    else
    {
        mpRefinedModelPart = mrOriginModelPart.CreateSubModelPart(refined_model_part_name);
    }
}

void NestedRefiningProcess::CheckRefiningCondition()
{
    auto gradient_variable = KratosComponents< Variable<array_1d<double, 3>>>::Get( mParameters["gradient_variable"].GetString());
    std::vector<double> scalar_values;
    ModelPart::NodesContainerType::iterator inodebegin = mpOwnModelPart->NodesBegin();
    for(IndexType i = 0; i < mpOwnModelPart->Nodes().size(); i++)
    {
        ModelPart::NodesContainerType::iterator inode = inodebegin+i;
        double value = std::log(norm_2(inode->FastGetSolutionStepValue(gradient_variable)));
        scalar_values.push_back(value);
    }
    
    double sum = std::accumulate(scalar_values.begin(), scalar_values.end(), 0.0);
    double mean = sum / scalar_values.size();

    std::vector<double> diff(scalar_values.size());
    std::transform(scalar_values.begin(), scalar_values.end(), diff.begin(), [mean](double x){return x -mean;});
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / scalar_values.size());

    std::lognormal_distribution<> d(mean, stdev);

    for(IndexType i = 0; i< mpOwnModelPart->Nodes().size(); i++)
    {
        ModelPart::NodesContainerType::iterator inode = inodebegin+i;
        //double p = 
    }
}

}  // namespace Kratos.
