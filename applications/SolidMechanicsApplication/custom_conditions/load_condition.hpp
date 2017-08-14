//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/boundary_condition.hpp"

#include "custom_utilities/solid_mechanics_math_utilities.hpp"

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

/// Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a General Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) LoadCondition
    : public BoundaryCondition
{
public:

    ///@name Type Definitions
    ///@{
  
    // Counted pointer of LoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( LoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    LoadCondition();
  
    /// Default constructor.
    LoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    LoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    LoadCondition( LoadCondition const& rOther);

    /// Destructor
    virtual ~LoadCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const override;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const override;


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

    /**
     * Initialize General Variables
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculate the External Load of the Condition
     */
    virtual void CalculateExternalLoad(GeneralVariables& rVariables);


    /**
     * Calculation of the External Forces Vector for a force or pressure vector 
     */
    virtual void CalculateAndAddExternalForces(Vector& rRightHandSideVector,
					       GeneralVariables& rVariables,
					       double& rIntegrationWeight ) override;


    /**
     * Calculation of the External Forces Vector for a force or pressure vector 
     */
    virtual double& CalculateAndAddExternalEnergy(double& rEnergy,
						  GeneralVariables& rVariables,
						  double& rIntegrationWeight ) override;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;


}; // class LoadCondition.

} // namespace Kratos.

#endif // KRATOS_LOAD_CONDITION_H_INCLUDED defined 
