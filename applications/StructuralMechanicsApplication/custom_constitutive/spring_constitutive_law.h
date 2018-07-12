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
template<std::size_t TDim>
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

    /**
     * Counted pointer of SpringConstitutiveLaw
     */
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
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_X );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_Y );
    KRATOS_DEFINE_LOCAL_FLAG( NULL_ROTATTIONAL_STIFFNESS_Z );

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
    ~SpringConstitutiveLaw() override;

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
    Pointer Create(Kratos::Parameters NewParameters) override
    {
        return Kratos::make_shared<SpringConstitutiveLaw<TDim>>(NewParameters);
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
    
    Flags mConstitutiveLawFlags; /// Constitutive flags

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
        rSerializer.save("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.save("TimeInterval",mTimeInterval);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLawFlags",mConstitutiveLawFlags);
        rSerializer.load("TimeInterval",mTimeInterval);
    }


}; // Class SpringConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_SPRING_LAW_H_INCLUDED  defined
