// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


// System includes

// External includes

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/spring_constitutive_law.h"

namespace Kratos
{
namespace Python
{

using namespace pybind11;


void  AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    class_< SpringConstitutiveLaw, typename SpringConstitutiveLaw::Pointer, ConstitutiveLaw >
    (m, "SpringConstitutiveLaw").def(init<>() )
    ;
}

}  // namespace Python.
}  // namespace Kratos.
