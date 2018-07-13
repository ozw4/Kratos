//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include "includes/define.h"
#include "utilities/math_utils.h"
#include "geometries/geometry.h"

// External includes

// Project includes
#include "custom_elements/iga_truss_element_3DXN.h"

// Application includes
#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

namespace Kratos {

IgaTrussElement3DXN::IgaTrussElement3DXN(
    IndexType NewId,
    GeometryType::Pointer pGeometry
)
    : CurveBaseDiscreteElement(NewId, pGeometry)
{
}

IgaTrussElement3DXN::IgaTrussElement3DXN(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
)
    : CurveBaseDiscreteElement(NewId, pGeometry, pProperties)
{
}

IgaTrussElement3DXN::IgaTrussElement3DXN(
)
    : CurveBaseDiscreteElement()
{
}

IgaTrussElement3DXN::~IgaTrussElement3DXN(
)
{
}

Element::Pointer
IgaTrussElement3DXN::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
) const
{
    return Kratos::make_shared<IgaTrussElement3DXN>(NewId,
        GetGeometry().Create(ThisNodes), pProperties);
};

void
IgaTrussElement3DXN::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    size_t numberOfNodes = GetGeometry().size();
    size_t numberOfDofs = numberOfNodes * 3;

    if (rElementalDofList.size() != numberOfDofs) {
        rElementalDofList.resize(numberOfDofs);
    }

    for (size_t i = 0; i < numberOfNodes; i++) {
        auto& node = GetGeometry()[i];

        size_t index = i * 3;

        rElementalDofList[index + 0] = node.pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = node.pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = node.pGetDof(DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
};

void
IgaTrussElement3DXN::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY;

    size_t numberOfNodes = GetGeometry().size();
    size_t numberOfDofs = numberOfNodes * 3;

    if (rResult.size() != numberOfDofs) {
        rResult.resize(numberOfDofs);
    }

    for (size_t i = 0; i < numberOfNodes; i++) {
        auto& node = GetGeometry()[i];

        size_t index = i * 3;

        rResult[index + 0] = node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = node.GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH("")
};

void
IgaTrussElement3DXN::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
)
{
}

void
IgaTrussElement3DXN::InitializeMaterial(
)
{
    KRATOS_TRY

    KRATOS_CATCH("");
}

void
IgaTrussElement3DXN::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
)
{
    KRATOS_TRY

    size_t numberOfNodes = GetGeometry().size();
    size_t numberOfDofs = numberOfNodes * 3;

    if (CalculateStiffnessMatrixFlag) {
        if (rLeftHandSideMatrix.size1() != numberOfDofs) {
            rLeftHandSideMatrix.resize(numberOfDofs, numberOfDofs);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix(numberOfDofs, numberOfDofs);
    }

    if (CalculateResidualVectorFlag) {
        if (rRightHandSideVector.size() != numberOfDofs) {
            rRightHandSideVector.resize(numberOfDofs);
        }

        noalias(rRightHandSideVector) = ZeroVector(numberOfDofs);
    }

    // get integration data
    
    double& integrationWeight = GetValue(INTEGRATION_WEIGHT);
    Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // get properties

    auto& properties = GetProperties();

    double E = properties[YOUNG_MODULUS];
    double A = properties[CROSS_AREA];
    double prestress = properties[PRESTRESS_CAUCHY];

    // compute base vectors

    Vector3d actualBaseVector = ZeroVector(3);
    GetBaseVector(actualBaseVector, DN_De);

    double referenceA = norm_2(mBaseVector0);
    double actualA = norm_2(actualBaseVector);

    // green-lagrange strain
    double E11_membrane = 0.5 * (actualA * actualA - referenceA * referenceA);

    // normal force
    double S11_membrane = prestress * A + E11_membrane * A * E / (referenceA *
        referenceA);

    // 1st variation of the axial strain 
    Vector epsilonVarDof = ZeroVector(numberOfDofs);
    Get1stVariationsAxialStrain(epsilonVarDof, actualBaseVector, 3, DN_De);
    epsilonVarDof = epsilonVarDof / (referenceA * referenceA);

    // 2nd variation of the axial strain 
    Matrix epsilonVar2Dof = ZeroMatrix(numberOfDofs, numberOfDofs);
    Get2ndVariationsAxialStrain(epsilonVar2Dof, 3, DN_De);
    epsilonVar2Dof = epsilonVar2Dof / (referenceA * referenceA);

    for (size_t r = 0; r < numberOfDofs; r++) {
        for (size_t s = 0; s < numberOfDofs; s++) {
            rLeftHandSideMatrix(r, s) = E * A * epsilonVarDof[r] *
                epsilonVarDof[s] + S11_membrane * epsilonVar2Dof(r, s);
        }
    }

    rRightHandSideVector = -S11_membrane * epsilonVarDof;

    rLeftHandSideMatrix *= integrationWeight;
    rRightHandSideVector *= integrationWeight;

    KRATOS_CATCH("");
}

void
IgaTrussElement3DXN::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if (rOutput.size() != 1) {
        rOutput.resize(1);
    }

    rOutput[0] = 0.0;
}

void
IgaTrussElement3DXN::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
)
{
    if (rValues.size() != 1) {
        rValues.resize(1);
    }

    rValues[0] = ZeroVector(3);
}

void
IgaTrussElement3DXN::save(
    Serializer& rSerializer
) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

void
IgaTrussElement3DXN::load(
    Serializer& rSerializer
)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

} // namespace Kratos
