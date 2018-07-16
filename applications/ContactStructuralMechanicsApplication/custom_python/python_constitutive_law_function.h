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

#if !defined (KRATOS_PYTHON_CONSTITUTIVE_LAW_FUNCTION)
#define  KRATOS_PYTHON_CONSTITUTIVE_LAW_FUNCTION

// System includes

// External includes

// Project includes
#include "includes/node.h"

namespace Kratos
{
namespace Python
{
///@name Python Globals
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
///@name Python Classes
///@{

/**
 * @class ConstitutiveLawFunction
 * @ingroup StructuralMechanicsApplication
 * @brief This function allows to call a function method of type f(x, y, z, t) implemented in python.
 * @details Uses python functions to evaluate the bahaviour of the CL
 * The functions can be constructed by providing a python-defined method of the type
 *
 *  class aux_object_cpp_callback:
 *    def __init__(self, function_string ):
 *        #TODO: check python version
 *        self.compiled_function = compile(function_string, '', 'eval', optimize=2)
 *
 *    def f(self,x,y,z,t):
 *        return  eval(self.compiled_function)
 *
 * the object is then insantiated as
 * aux_function = ConstitutiveLawFunction(aux_object_cpp_callback(self.function_string))
 * @note Based on python_function_callback_utility
 * @note This makes this file to depend on python
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class ConstitutiveLawFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// The node defintion
    typedef Node<3> NodeType;

    /// Counted pointer of ConstitutiveLawFunction
    KRATOS_CLASS_POINTER_DEFINITION(ConstitutiveLawFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The text containing the function
     */
    ConstitutiveLawFunction( const std::string& rFunctionBody);

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This methods returns the value of the function
     * @param ThisNode The node whre to evaluate
     * @param Time The current time
     * @return The resulting value of the function
     */
    double CallFunction(
        const NodeType& ThisNode,
        std::vector<std::string>& AdditionalVariables,
        const double Time
        );

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

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
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; // ConstitutiveLawFunction
}  // namespace Python.
}  // namespace Kratos.
#endif // KRATOS_PYTHON_CONSTITUTIVE_LAW_FUNCTION  defined
