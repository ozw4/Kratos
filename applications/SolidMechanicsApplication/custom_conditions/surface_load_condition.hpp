//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/load_condition.hpp"

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

// Surface load condition for 3D geometries.

class KRATOS_API(SOLID_MECHANICS_APPLICATION) SurfaceLoadCondition
    : public LoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of SurfaceLoadCondition
    KRATOS_CLASS_POINTER_DEFINITION( SurfaceLoadCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SurfaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    SurfaceLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    SurfaceLoadCondition( SurfaceLoadCondition const& rOther);

    /// Destructor
    virtual ~SurfaceLoadCondition();

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
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;



    //************* COMPUTING  METHODS


    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

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
    SurfaceLoadCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize System Matrices
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculation of the Position Increment
     */
    virtual Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Calculation of the Total Position Increment
     */
    virtual Matrix& CalculateTotalDeltaPosition(Matrix & rDeltaPosition);


    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, 
				     const double& rPointNumber);

    /**
     * Calculate the External Load of the Condition
     */
    virtual void CalculateExternalLoad(GeneralVariables& rVariables);

    /**
     * Calculation of the Load Stiffness Matrix which usually is subtracted to the global stiffness matrix
     */
    virtual void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				     GeneralVariables& rVariables,
				     double& rIntegrationWeight);


    //utilities::

    void MakeCrossMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
			 Vector& U );

    void CrossProduct(Vector& cross,
		      Vector& a,
		      Vector& b );


    void AddMatrix(MatrixType& Destination,
		   boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
		   int InitialRow,
		   int InitialCol );

    void SubtractMatrix(MatrixType& Destination,
			boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
			int InitialRow,
			int InitialCol );


    void ExpandReducedMatrix(Matrix& Destination,
			     Matrix& ReducedMatrix );


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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

}; // class SurfaceLoadCondition.

} // namespace Kratos.

#endif // KRATOS_SURFACE_LOAD_CONDITION_H_INCLUDED defined 
