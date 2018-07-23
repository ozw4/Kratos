//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//
// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"


namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

}

}  // namespace Python.
} // Namespace Kratos
