// =============================================================================================
// postprocess.cc  (UPDATED for NEW ICE APP: u, phi ONLY)
// =============================================================================================
//
// Changes
//
// original postprocessing file was written for the old CH–AC model with
// primary fields:
//    variable 0 = c
//    variable 1 = n
//
// In the ice app the primary fields are:
//    variable 0 = u
//    variable 1 = phi
//
// constants (c, n, Kn, etc.).  must be updated or removed
//
// INCREMENTAL STRATEGY 
///   * only depend on existing variables (u, phi)
///   * do not require new Demange constants yet
//
// NEW POSTPROCESS OUTPUTS :
//   Variable 0: "phi_interface"
//   Variable 1: "u_value"



// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{
  // ---------------------------------------------------------------------------
  // Postprocess Variable 0: phi_interface = (1 - phi^2) interface indicator
  // ---------------------------------------------------------------------------
  set_variable_name(0, "phi_interface");
  set_variable_type(0, SCALAR);

  // Depends only on phi value 
  set_dependencies_value_term_RHS(0, "phi");
  set_dependencies_gradient_term_RHS(0, "");

  set_output_integral(0, true);

  // ---------------------------------------------------------------------------
  // Postprocess Variable 1: u_value = u
  // ---------------------------------------------------------------------------
  set_variable_name(1, "u_value");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "u");
  set_dependencies_gradient_term_RHS(1, "");

  set_output_integral(1, false);
}

//original*

// =============================================================================================
// postProcessedFields: Set the postprocessing expressions
// =============================================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::postProcessedFields(
  [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>>
    &variable_list,
  [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>
                                                            &pp_variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<double>> q_point_loc,
  [[maybe_unused]] const VectorizedArray<double>             element_volume) const
{
  // ---------------------------------------------------------------------------
  // Get primary field values 
  //   variable 0 = u
  //   variable 1 = phi
  // ---------------------------------------------------------------------------

  scalarvalueType u   = variable_list.get_scalar_value(0);
  scalarvalueType phi = variable_list.get_scalar_value(1);

  // ---------------------------------------------------------------------------
  // Compute postprocess expressions 
  // ---------------------------------------------------------------------------

  //     Interface indicator:
  //     phi_interface = 1 - phi^2
  scalarvalueType phi_interface = constV(1.0) - phi * phi;

  // u itself
  scalarvalueType u_value = u;

  // ---------------------------------------------------------------------------
  // Submitting the terms for the postprocessing expressions (no gradient term)
  // ---------------------------------------------------------------------------
  pp_variable_list.set_scalar_value_term_RHS(0, phi_interface);
  pp_variable_list.set_scalar_value_term_RHS(1, u_value);
}
