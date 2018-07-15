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

#if !defined (KRATOS_SPRING_LAW_H_INCLUDED)
#define  KRATOS_SPRING_LAW_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
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

/**
 * @class PythonConstitutiveLawFunction
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
 * aux_function = PythonConstitutiveLawFunction(aux_object_cpp_callback(self.function_string))
 * @note Based on python_function_callback_utility
 * @note This makes this file to depend on python
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class PythonConstitutiveLawFunction
{
public:
    ///@name Type Definitions
    ///@{

    /// The node defintion
    typedef Node<3> NodeType;

    /// Counted pointer of PythonConstitutiveLawFunction
    KRATOS_CLASS_POINTER_DEFINITION(PythonConstitutiveLawFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rFunctionBody The text containing the function
     */
    PythonConstitutiveLawFunction( const std::string& rFunctionBody);

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
}; // PythonConstitutiveLawFunction

/**
 * @class SpringConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief Spring/damper/mass/inertia... constitutive law for 3D and 2D points
 * @details Uses python functions to evaluate the bahaviour of the CL
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SpringConstitutiveLaw
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The process info type
    typedef ProcessInfo       ProcessInfoType;

    /// The constitutive law definiton
    typedef ConstitutiveLaw          BaseType;

    /// The size definition
    typedef std::size_t              SizeType;

    /// The index type
    typedef std::size_t             IndexType;

    /// Counted pointer of SpringConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION( SpringConstitutiveLaw );

    /**
     * @brief Flags related to constitutive law computation
     */
    KRATOS_DEFINE_LOCAL_FLAG( NULL_MASS );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_INERTIA_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_STIFFNESS_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_STIFFNESS_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_DAMPING_RATIO_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_DAMPING_RATIO_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_DAMPING_RATIO_Z );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_DAMPING_RATIO_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_DAMPING_RATIO_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATIONAL_DAMPING_RATIO_Z );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_MASS );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_INERTIA_X );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_INERTIA_Y );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_INERTIA_Z );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_STIFFNESS_Z );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_STIFFNESS_Z );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_DAMPING_RATIO_X );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_DAMPING_RATIO_Y );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_DAMPING_RATIO_Z );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_DAMPING_RATIO_X );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_DAMPING_RATIO_Y );
    KRATOS_DEFINE_LOCAL_FLAG( CONSTANT_ROTATIONAL_DAMPING_RATIO_Z );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    SpringConstitutiveLaw();

    /**
     * @brief Default constructor. Using parameters instead of the one by default
     * @param NewParameters The configuration parameters of the new constitutive law
     */
    SpringConstitutiveLaw(Kratos::Parameters NewParameters);

    /**
     * @brief Copy constructor.
     */
    SpringConstitutiveLaw (const SpringConstitutiveLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~SpringConstitutiveLaw() override {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * @note implementation scheme:
     *      ConstitutiveLaw::Pointer p_clone(new ConstitutiveLaw());
     *      return p_clone;
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override
    {
        return Kratos::make_shared<SpringConstitutiveLaw>(NewParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 1;
    }
    
    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    array_1d<double, 3 > & CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

protected:

    ///@name Protected static Member Variables
    ///@{
    
    ///@}
    ///@name Protected member Variables
    ///@{
    
    /* The definition of the functions */
    std::unordered_map<IndexType, PythonConstitutiveLawFunction::Pointer> mFunctions; /// This map stores the defined functions
    std::unordered_map<IndexType, double> mConstantValues;                            /// This method stores the constant values

    std::vector<std::string> mAdditionalDependenceVariables;                          /// Variables to include additional dependence

    IndexType mNodalIndex = 0;          /// The index of the current node on the geometry

    Flags mConstitutiveLawFlags;       /// Constitutive flags

    array_1d<double, 2> mTimeInterval; /// The time interval

    ///@}
    ///@name Protected Operators
    ///@{
    
    ///@}
    ///@name Protected Operations
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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
//         rSerializer.save("Functions",mFunctions);
        rSerializer.save("ConstantValues",mConstantValues);
        rSerializer.save("AdditionalDependenceVariables",mAdditionalDependenceVariables);
        rSerializer.save("NodeIndex",mNodalIndex);
        rSerializer.save("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.save("TimeInterval",mTimeInterval);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
//         rSerializer.load("Functions",mFunctions);
        rSerializer.load("ConstantValues",mConstantValues);
        rSerializer.load("AdditionalDependenceVariables",mAdditionalDependenceVariables);
        rSerializer.load("NodeIndex",mNodalIndex);
        rSerializer.load("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.load("TimeInterval",mTimeInterval);
    }


}; // Class SpringConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_SPRING_LAW_H_INCLUDED  defined
