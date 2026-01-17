// ===========================================================================
// icsbcs.cc (ICE MODEL SCAFFOLD: u, phi ONLY)
// ===========================================================================
//
// changes in this version
//
// - REMOVED c/n initial conditions
//
// - ADDED simple ICs:
//     u   = 0 everywhere (safe default)
//     phi = -1 , +1 in seed(s) using tanh profile (smooth interface)
//
// - BC function left blank 
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
  // same interface width style as original code
  //removed concentration values
  const double interface_width = 1.0;
  //Default value
  scalar_IC = 0.0;

  // -----------------------------
  // u IC (variable 0)
  // -----------------------------
  if (index == 0)
    {
      // Minimal safe default for incremental runs
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
      //Computes distance from the current point p to center1
      //Euclidean distance formula
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

      //combine seeds
      // Use max(s1,s2) to avoid phi exceeding +1 if seeds overlap
      double s = s1;
      if (s2 > s)
        s = s2;

      // Map s in [0,1] -> phi in [-1,+1]
      //s=1 → phi=+1 (inside a seed)
      //s=0 → phi=-1 (outside)
      scalar_IC = 2.0 * s - 1.0;
      return;
    }

  // Any unexpected index (shouldn't happen if num vars = 2)
  scalar_IC = 0.0;
}

////original*

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
  //
  //left blank for incremental step.
  //
}
////original*
