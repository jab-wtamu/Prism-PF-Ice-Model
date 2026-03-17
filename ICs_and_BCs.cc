// ===========================================================================
// icsbcs.cc (ICE MODEL SCAFFOLD: u, phi ONLY)
// ===========================================================================
  // VISUAL VERIFICATION IC FOR u: Gaussian “bump” initial condition
  // ===========================================================================
  
  //
  // MATHEMATICAL FORM: Radial Gaussian
  //   u(x) = u0 + amp * exp( -|x - x0|^2 / (2*sigma^2) )
  //
  // PARAMETERS :
  // - u0    : baseline value of u far from the bump
  // - amp   : bump height (bigger -> stronger gradients -> faster visible change)
  // - sigma : bump width (smaller -> sharper bump -> stronger initial gradients)
  
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
  // same interface width style as  original code
  const double interface_width = 1.0;

  scalar_IC = 0.0;

  // -----------------------------
  // u IC (variable 0)
  // -----------------------------
   if (index == 0)
{
  const bool use_gaussian_bump = false;  // <-- flip this when ready

const double base = u_initial;  // from parameters.prm

  if (!use_gaussian_bump)
  {
    scalar_IC = base;   // uniform supersaturation
    return;
  }

  const double amp   = 0.2;
  const double sigma = 8.0;

  double r2 = 0.0;
  for (unsigned int dir = 0; dir < dim; ++dir)
  {
    const double dx = p[dir] - center1[dir];
    r2 += dx * dx;
  }

scalar_IC = base + amp * std::exp(-r2 / (2.0 * sigma * sigma));
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
scalar_IC = phi_matrix + (phi_seed - phi_matrix) * s;      return;
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
  if (index == 0) // u
{
  scalar_BC = u0;   // far-field supersaturation from parameters.prm
  return;
}

}
