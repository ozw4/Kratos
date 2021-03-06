//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#if defined(KRATOS_PYTHON)

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_python/add_cross_sections_to_python.h"

#include "solid_mechanics_application.h"

namespace Kratos
{

namespace Python
{

using namespace pybind11;


PYBIND11_MODULE(KratosSolidMechanicsApplication,m)
{

  class_<KratosSolidMechanicsApplication,
         KratosSolidMechanicsApplication::Pointer,
         KratosApplication>(m,"KratosSolidMechanicsApplication")
      .def(init<>())
      ;

  AddCustomUtilitiesToPython(m);
  AddCustomStrategiesToPython(m);
  AddCustomConstitutiveLawsToPython(m);
  AddCrossSectionsToPython(m);
  AddCustomProcessesToPython(m);

  //registering variables in python ( if must to be seen from python )

  // Generalized eigenvalue problem
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BUILD_LEVEL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EIGENVALUE_VECTOR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EIGENVECTOR_MATRIX )

  // For process information
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VECTOR_TIME_INTEGRATION_METHODS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPONENT_TIME_INTEGRATION_METHODS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_TIME_INTEGRATION_METHODS )

  // For integration schemes
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COMPUTE_CONSISTENT_MASS_MATRIX )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, STEP_DISPLACEMENT )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, STEP_ROTATION )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MIDDLE_VELOCITY )

  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, POINT_LOAD )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, FORCE_LOAD )

  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, MOMENT_LOAD )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ELASTIC_LOAD )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FORCE_LOAD_VECTOR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, MOMENT_LOAD_VECTOR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ELASTIC_LOAD_VECTOR )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POSITIVE_FACE_PRESSURE_VECTOR )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NEGATIVE_FACE_PRESSURE_VECTOR )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLANE_MOMENT_LOAD )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PLANE_MOMENT_LOAD_VECTOR )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BALLAST_COEFFICIENT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, BALLAST_COEFFICIENT_VECTOR )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, WRITE_ID )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TIME_INTEGRATION_ORDER )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_ALPHA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_BETA )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DISPLACEMENT_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROTATION_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_REACTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VON_MISES_STRESS )

  // For beams
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_SECTION_AREA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_SECTION_RADIUS )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_SECTION_SIDES )

  // For shells cross section
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_CROSS_SECTION )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

  // For shell generalized variables
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_STRAIN )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_STRAIN_GLOBAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_CURVATURE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_CURVATURE_GLOBAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_FORCE )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_FORCE_GLOBAL )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_MOMENT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHELL_MOMENT_GLOBAL )

  // Reading beam section properties
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SECTION_HEIGHT )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SECTION_WIDTH  )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INERTIA_X )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, INERTIA_Y )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SECTION_SIZE )

  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNGxAREA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNGxINERTIA_X )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, YOUNGxINERTIA_Y )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEARxREDUCED_AREA )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEARxPOLAR_INERTIA )

  }


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
