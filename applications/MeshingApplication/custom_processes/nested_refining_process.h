//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		     BSD License
//					         Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined( KRATOS_NESTED_REFINING_PROCESS_H_INCLUDED )
#define  KRATOS_NESTED_REFINING_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>


// External includes


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"


namespace Kratos {

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

/// Generate a nested mesh in the region of interest
/**
 * The process creates a sub model part with the nested mesh
 * This process can be called within the refined sub model part
 * The origin model part is stored under a new sub model part (own):
 * MainModelPart
 *     Own
 *     RefinedSubModelPart
 *         Own
 *         RefinedSubModelPart
 *             Own
 *             ...
 */
class NestedRefiningProcess : public Process {
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NestedRefiningProcess
    KRATOS_CLASS_POINTER_DEFINITION(NestedRefiningProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    NestedRefiningProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
    );

    /// Destructor.
    virtual ~NestedRefiningProcess() {
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

    void Execute() override;

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
        return "NestedRefiningProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NestedRefiningProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ModelPart& mrModelPart;
    Parameters mParameters;

    unsigned int mEchoLevel;

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
    NestedRefiningProcess& operator=(NestedRefiningProcess const& rOther);

    /// Copy constructor.
    //NestedRefiningProcess(NestedRefiningProcess const& rOther);

    ///@}

}; // Class NestedRefiningProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  NestedRefiningProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const NestedRefiningProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_NESTED_REFINING_PROCESS_H_INCLUDED
