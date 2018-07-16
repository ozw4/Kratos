import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignNodalElementsToNodes(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignNodalElementsToNodes(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                        : 0,
            "main_model_part"                : "Structure",
            "sub_model_part_name"            : "",
            "rayleigh_damping"               : false,
            "assign_active_flag_node"        : true,
            "nodal_mass"                     : null,
            "nodal_inertia"                  : [null, null, null],
            "nodal_stiffness"                : [null, null, null],
            "nodal_rotational_stiffness"     : [null, null, null],
            "nodal_damping_ratio"            : [null, null, null],
            "nodal_rotational_damping_ratio" : [null, null, null],
            "additional_dependence_variables": [],
            "interval"                       : [0.0, 1e30]
        }
        """
        )

        # Overwrite the default settings with user-provided parameters
        self.settings = settings
        self.settings.RecursivelyValidateAndAssignDefaults(default_settings)

        # The main model part
        self.model = Model
        self.main_model_part = self.model[self.settings["model_part_name"].GetString()]

        # The creation of the process
        process_parameters = KratosMultiphysics.Parameters("""{}""")
        process_parameters.AddValue("main_model_part", self.settings["sub_model_part_name"])
        process_parameters.AddValue("rayleigh_damping", self.settings["rayleigh_damping"])
        process_parameters.AddValue("assign_active_flag_node", self.settings["assign_active_flag_node"])
        process_parameters.AddValue("nodal_mass", self.settings["nodal_mass"])
        process_parameters.AddValue("nodal_inertia", self.settings["nodal_inertia"])
        process_parameters.AddValue("nodal_stiffness", self.settings["nodal_stiffness"])
        process_parameters.AddValue("nodal_rotational_stiffness", self.settings["nodal_rotational_stiffness"])
        process_parameters.AddValue("nodal_damping_ratio", self.settings["nodal_damping_ratio"])
        process_parameters.AddValue("nodal_rotational_damping_ratio", self.settings["nodal_rotational_damping_ratio"])
        process_parameters.AddValue("additional_dependence_variables", self.settings["additional_dependence_variables"])
        process_parameters.AddValue("interval", self.settings["interval"])
        self.assign_nodal_elements_to_nodes = StructuralMechanicsApplication.AssignNodalElementsToNodes(self.main_model_part, process_parameters)

    def ExecuteInitialize(self):
        self.assign_nodal_elements_to_nodes.ExecuteInitialize()
