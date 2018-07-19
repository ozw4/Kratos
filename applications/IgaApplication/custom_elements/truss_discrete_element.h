//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED)
#define KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/element.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

// Project includes
#include "custom_elements/curve_base_discrete_element.h"

namespace Kratos {

/**
 * @class TrussDiscreteElement
 *
 * @brief This is a 3D-X-node isogeometric truss element with 3 translational
 * dofs per node
 */
class TrussDiscreteElement
    : public CurveBaseDiscreteElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(TrussDiscreteElement);


    TrussDiscreteElement(
    );

    TrussDiscreteElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry
    );

    TrussDiscreteElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
    );

    ~TrussDiscreteElement() override;


    Element::Pointer
    Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
    * @brief This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void
    CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) override;

    /**
    * @brief Calculate a double Variable on the Element Constitutive Law
    * @param rVariable: The variable we want to get
    * @param rOutput: The values obtained int the integration points
    * @param rCurrentProcessInfo: the current process info instance
    */
    void
    CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Calculate a Vector Variable on the Element Constitutive Law
    * @param rVariable: The variable we want to get
    * @param rOutput: The values obtained int the integration points
    * @param rCurrentProcessInfo: the current process info instance
    */
    void
    CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rResult the ID's of the element degrees of freedom
    * @param rResult: The vector containing the equation id
    * @param rCurrentProcessInfo: The current process info instance
    */
    void
    EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
    * @param rElementalDofList: The vector containing the dof of the element
    * @param rCurrentProcessInfo: The current process info instance
    */
    void
    GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    ) override;

protected:
    /**
    * It initializes the material
    */
    virtual void
    InitializeMaterial(
    );

    /**
    * Called at the end of eahc solution step
    * @param rCurrentProcessInfo: the current process info instance
    */
    void
    FinalizeSolutionStep(
        ProcessInfo& rCurrentProcessInfo
    ) override;

private:
    friend class Serializer;

    void
    save(
        Serializer& rSerializer
    ) const override;

    void
    load(
        Serializer& rSerializer
    ) override;
}; // class TrussDiscreteElement

} // namespace Kratos

#endif // !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED)
