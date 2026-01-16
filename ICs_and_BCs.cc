// ===========================================================================
// icsbcs.cc (ICE MODEL SCAFFOLD: u, phi ONLY)
// ===========================================================================
//
// WHAT CHANGED (and WHY)
//
// - REMOVED c/n initial condition branching
//   WHY: new ice app has only u and phi
//
// - ADDED simple ICs:
//     u   = 0 everywhere (safe default; replace later with your u_infty)
//     phi = -1 in matrix, +1 in seed(s) using tanh profile (smooth interface)
//
// - BC function left blank to preserve "no non-uniform Dirichlet" behavior.
//   If your parameter file requests non-uniform Dirichlet BCs for u or phi,
//   we must implement them here.
//
// ===========================================================================


// ===========================================================================
// FUNCTION FOR INITIAL CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                                            [[maybe_unused]] const unsigned int index,
                                            [[maybe_unused]] double            &scalar_IC,
                                            [[maybe_unused]] Vector<double>    &vector_IC)
{
  // same interface width style as your original code
  const double interface_width = 1.0;

  scalar_IC = 0.0;

  // -----------------------------
  // u IC (variable 0)
  // -----------------------------
  if (index == 0)
    {
      // Minimal safe default for incremental runs
      // Later: set to far-field supersaturation (e.g., u_infty)
      scalar_IC = 0.0;
      return;
    }

  // -----------------------------
  // phi IC (variable 1)
  // -----------------------------
  if (index == 1)
    {
      // Seed 1
      double dist1 = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        dist1 += (p[dir] - center1[dir]) * (p[dir] - center1[dir]);
      dist1 = std::sqrt(dist1);

      const double s1 = 0.5 * (1.0 - std::tanh((dist1 - radius1) / interface_width));

      // Seed 2
      double dist2 = 0.0;
      for (unsigned int dir = 0; dir < dim; dir++)
        dist2 += (p[dir] - center2[dir]) * (p[dir] - center2[dir]);
      dist2 = std::sqrt(dist2);

      const double s2 = 0.5 * (1.0 - std::tanh((dist2 - radius2) / interface_width));

      // Use max(s1,s2) to avoid phi exceeding +1 if seeds overlap
      double s = s1;
      if (s2 > s)
        s = s2;

      // Map s in [0,1] -> phi in [-1,+1]
      scalar_IC = 2.0 * s - 1.0;
      return;
    }

  // Any unexpected index (shouldn't happen if num vars = 2)
  scalar_IC = 0.0;
}


// ===========================================================================
// FUNCTION FOR NON-UNIFORM DIRICHLET BOUNDARY CONDITIONS
// ===========================================================================

template <int dim, int degree>
void
customPDE<dim, degree>::setNonUniformDirichletBCs(
  [[maybe_unused]] const Point<dim>  &p,
  [[maybe_unused]] const unsigned int index,
  [[maybe_unused]] const unsigned int direction,
  [[maybe_unused]] const double       time,
  [[maybe_unused]] double            &scalar_BC,
  [[maybe_unused]] Vector<double>    &vector_BC)
{
  // Intentionally left blank for incremental step.
  //
  // If later you want (example) u = u_infty on all outer boundaries:
  //   if (index == 0) scalar_BC = u_infty;
  //
  // Or if you want phi fixed on boundaries, implement index == 1 here.
}
