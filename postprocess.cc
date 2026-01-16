// =============================================================================================
// postprocess.cc  (UPDATED for NEW ICE APP: u, phi ONLY)
// =============================================================================================
//
// WHAT CHANGED (and WHY)
//
// Your original postprocessing file was written for the old CH–AC model with
// primary fields:
//    variable 0 = c
//    variable 1 = n
//
// It computed:
//   (1) f_tot  : total free energy density (chemical + gradient)
//   (2) c_grad : |grad(c)|
//
// In the NEW ice app the primary fields are:
//    variable 0 = u
//    variable 1 = phi
//
// Therefore the old postprocess expressions referenced non-existent fields and
// constants (c, n, Kn, etc.). This must be updated or removed, otherwise the
// app will fail to compile/run when postprocessing is enabled.
//
// INCREMENTAL STRATEGY (to keep it running):
// - Provide simple, safe postprocessing quantities that:
///   * only depend on existing variables (u, phi)
///   * do not require new Demange constants yet
// - Keep the same "two postprocess vars" structure to minimize framework changes.
//
// NEW POSTPROCESS OUTPUTS (minimal + useful):
//   Variable 0: "phi_interface"
//     - A scalar that highlights the interface region.
//     - We output (1 - phi^2), which is ~0 in bulk (phi≈±1) and peaks near phi=0.
//     - Useful for tracking interface without gradients.
//
//   Variable 1: "u_value"
//     - Just u itself (handy for quick sanity checks).
//
// NOTE:
// - We set dependencies accordingly.
// - We keep integrals consistent: you can choose to integrate phi_interface
//   to track total interface measure; u_value integral usually not needed.
//
// =============================================================================================


// =================================================================================
// Set the attributes of the postprocessing variables
// =================================================================================

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{
  // ---------------------------------------------------------------------------
  // Postprocess Variable 0: phi_interface = (1 - phi^2)
  // ---------------------------------------------------------------------------
  set_variable_name(0, "phi_interface");
  set_variable_type(0, SCALAR);

  // Depends only on phi value (no gradients needed for this minimal output)
  set_dependencies_value_term_RHS(0, "phi");
  set_dependencies_gradient_term_RHS(0, "");

  // Optional: integral gives a rough measure of total interface "amount"
  set_output_integral(0, true);

  // ---------------------------------------------------------------------------
  // Postprocess Variable 1: u_value = u
  // ---------------------------------------------------------------------------
  set_variable_name(1, "u_value");
  set_variable_type(1, SCALAR);

  set_dependencies_value_term_RHS(1, "u");
  set_dependencies_gradient_term_RHS(1, "");

  // Usually you do NOT need an integral of u, but keep false as a safe default
  set_output_integral(1, false);
}


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
  // Get primary field values (NEW ice app indexing)
  //   variable 0 = u
  //   variable 1 = phi
  // ---------------------------------------------------------------------------

  scalarvalueType u   = variable_list.get_scalar_value(0);
  scalarvalueType phi = variable_list.get_scalar_value(1);

  // ---------------------------------------------------------------------------
  // Compute postprocess expressions (minimal, no extra constants)
  // ---------------------------------------------------------------------------

  // (1) Interface indicator: ~0 in bulk, peaks at interface
  //     phi_interface = 1 - phi^2
  scalarvalueType phi_interface = constV(1.0) - phi * phi;

  // (2) u itself
  scalarvalueType u_value = u;

  // ---------------------------------------------------------------------------
  // Submit postprocessed scalars (value-term only; no gradient term)
  // ---------------------------------------------------------------------------
  pp_variable_list.set_scalar_value_term_RHS(0, phi_interface);
  pp_variable_list.set_scalar_value_term_RHS(1, u_value);
}
