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
#include "custom_processes/multi_scale_refining_process.h"
#include "utilities/sub_model_parts_list_utility.h"
#include "custom_utilities/uniform_refine_utility.h"

namespace Kratos
{

MultiScaleRefiningProcess::MultiScaleRefiningProcess(
    Model& rThisModel,
    Parameters ThisParameters
) : mrModel(rThisModel)
  , mParameters(ThisParameters)
{
    Parameters DefaultParameters = Parameters(R"(
    {
        "model_part_name"              : "MainModelPart",
        "echo_level"                   : 0,
        "maximum_number_of_levels"     : 4,
        "number_of_divisions_at_level"  : 2,
        "refining_boundary_condition"  : "Condition2D3N"
    }
    )");

    mParameters.ValidateAndAssignDefaults(DefaultParameters);

    mModelPartsNames.push_back(mParameters["model_part_name"].GetString());
    mEchoLevel = mParameters["echo_level"].GetInt();
    mMaxMultiScaleLevels = mParameters["maximum_number_of_levels"].GetInt();
    mDivisions = mParameters["number_of_divisions_at_level"].GetInt();
    mConditionName = mParameters["refining_boundary_condition"].GetString();

    InitializeModelParts();

}


void MultiScaleRefiningProcess::InitializeModelParts()
{
    ModelPart& root_mp = mrModel.GetModelPart(mModelPartsNames[0]);
    
    // StringVectorType sub_model_parts_names = root_mp.GetSubModelPartNames();
    StringVectorType sub_model_parts_names = this->GetRecursiveSubModelPartNames(root_mp);

    KRATOS_WATCH(mrModel)

    for (unsigned int i = 1; i <= mMaxMultiScaleLevels; i++)
    {
        std::string name = mModelPartsNames[0] + "-level_" + std::to_string(i);
        KRATOS_INFO_IF(this->Info(), mEchoLevel>1) << "Creating Model Part " << name << std::endl;
        mModelPartsNames.push_back(name);

        ModelPart::Pointer new_model_part = Kratos::make_shared<ModelPart>( name );
        mrModel.AddModelPart(new_model_part);
        KRATOS_WATCH(mrModel)
        // TODO: using the new Model interface:
        // mrModel.CreateModelPart(name);

        // for (auto sub_full_name : sub_model_parts_names)
        // {
        //     ModelPart::Pointer sub_model_part = new_model_part;
        //     KRATOS_WATCH(sub_full_name)
        //     std::istringstream iss(sub_full_name);
        //     std::string token;
        //     while (std::getline(iss, token, '.'))
        //     {
        //         if (sub_model_part->HasSubModelPart(token))
        //             sub_model_part = sub_model_part->pGetSubModelPart(token);
        //         else
        //             sub_model_part = sub_model_part->CreateSubModelPart(token);
        //     }
        //     KRATOS_WATCH(*new_model_part)
        // }
    }

    KRATOS_WATCH(" ")
    KRATOS_WATCH("Finished constructing Model Parts")
    KRATOS_WATCH(" ")

    KRATOS_WATCH(mrModel)

    for (auto name : mModelPartsNames)
    {
        ModelPart& mp = mrModel.GetModelPart(name);
        KRATOS_WATCH(name)
        KRATOS_WATCH(mp)
    }

    // for (auto sub_name : sub_model_parts_names)
    // {
    //     KRATOS_WATCH(sub_name)
    // }

}


StringVectorType MultiScaleRefiningProcess::GetRecursiveSubModelPartNames(ModelPart& rThisModelPart, std::string Prefix)
{
    StringVectorType names = rThisModelPart.GetSubModelPartNames();
    if (!Prefix.empty())
        Prefix += ".";
    
    for (auto& name : names)
    {
        ModelPart& sub_model_part = rThisModelPart.GetSubModelPart(name);
        auto sub_names = this->GetRecursiveSubModelPartNames(sub_model_part, Prefix + name);
        name.insert(0, Prefix);
        for (auto sub_name : sub_names)
            names.push_back(sub_name);
    }

    return names;
}

}  // namespace Kratos.