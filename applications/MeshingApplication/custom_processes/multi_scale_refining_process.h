//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, @{KRATOS_APP_AUTHOR}
//

#if !defined(KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED )
#define  KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"


namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

/**
 * Vector of strings type
 */
typedef std::vector<std::string> StringVectorType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Generate a new mesh in the region of interest
/**
 * The process creates another model part with the refined mesh
 * This process can be called within the refined model part
 * The origin model part and the refined model parts are stored on the same model:
 * Model
 *     MainModelPart
 *         SubModelParts
 *         ...
 *     MainModelPart-level_1
 *         SubModelParts
 *         ...
 *     MainModelPart-level_2
 *         SubModelParts
 *         ...
 */
class MultiScaleRefiningProcess : public Process {
public:

    ///@name Type Definitions
    ///@{

    /**
     * Node type
     */
    typedef Node<3> NodeType;

    /**
     * Type of IDs
     */
    typedef std::size_t IndexType;

    /**
     * Map types to locate nodes in the mesh
     */
    typedef std::unordered_map<IndexType, NodeType::Pointer> IndexNodeMapType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MultiScaleRefiningProcess
    KRATOS_CLASS_POINTER_DEFINITION(MultiScaleRefiningProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    MultiScaleRefiningProcess(
        Model& rThisModel,
        Parameters ThisParameters = Parameters(R"({})")
    );

    /// Destructor.
    virtual ~MultiScaleRefiningProcess() {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()() {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() override
    {}

    ModelPart& GetModelPartAtLevel(int Level);

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
    std::string Info() const override
    {
        return "MultiScaleRefiningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MultiScaleRefiningProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {}

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

    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    Parameters mParameters;

    unsigned int mEchoLevel;
    unsigned int mMaxMultiScaleLevels;
    unsigned int mDivisions;
    unsigned int mWorkingLevel;

    std::string mConditionName;

    StringVectorType mModelPartsNames;  /// The names of the model parts at each level

    std::vector<IndexNodeMapType> mCoarseRefinedNodesMap; /// Mapping from own to refined
    std::vector<IndexNodeMapType> mRefinedCoarseNodesMap; /// Mapping from refined to own

    ///@}
    ///@name Private Operators
    ///@{

    void InitializeModelParts();

    void SetWorkingLevel(const int Level);

    void InterpolateLevelBoundaryValuesAtSubStep(const int SubStep, const int SubSteps);

    void UpdateSubLevel();

    void TransferDataToCoarseLevel();

    StringVectorType GetRecursiveSubModelPartNames(ModelPart& rThisModelPart, std::string Prefix = "");

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
    // MultiScaleRefiningProcess& operator=(MultiScaleRefiningProcess const& rOther);

    /// Copy constructor.
    // MultiScaleRefiningProcess(MultiScaleRefiningProcess const& rOther);

    ///@}

}; // Class MultiScaleRefiningProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MultiScaleRefiningProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MultiScaleRefiningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_MULTI_SCALE_REFINING_PROCESS_H_INCLUDED defined
