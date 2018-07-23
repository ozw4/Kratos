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
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/apply_component_table_process.hpp"
#include "custom_processes/apply_double_table_process.hpp"
#include "custom_processes/apply_constant_hydrostatic_pressure_process.hpp"
#include "custom_processes/apply_hydrostatic_pressure_table_process.hpp"


namespace Kratos
{

namespace Python
{

void  AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<ApplyComponentTableProcess, ApplyComponentTableProcess::Pointer, Process>
    (m, "ApplyComponentTableProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyDoubleTableProcess, ApplyDoubleTableProcess::Pointer, Process>
    (m, "ApplyDoubleTableProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyConstantHydrostaticPressureProcess, ApplyConstantHydrostaticPressureProcess::Pointer, Process>
    (m, "ApplyConstantHydrostaticPressureProcess")
    .def(init < ModelPart&, Parameters>());
    class_<ApplyHydrostaticPressureTableProcess, ApplyHydrostaticPressureTableProcess::Pointer, Process>
    (m, "ApplyHydrostaticPressureTableProcess")
    .def(init < ModelPart&, Parameters>());
}

}  // namespace Python.
} // Namespace Kratos
