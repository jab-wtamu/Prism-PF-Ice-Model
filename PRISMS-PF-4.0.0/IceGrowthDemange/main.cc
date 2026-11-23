// Snow crystal growth (Demange et al.) -- PRISMS-PF 4.0.0, 3D only.

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  // Replaces "set dim" and "set degree" from the PF3 parameters file.
  constexpr unsigned int dim    = 3;
  constexpr unsigned int degree = 2;

  /**
   * Three fields, with the same names and indices as the PF3 application:
   *   u   (0) - reduced supersaturation
   *   phi (1) - ice/vapor order parameter
   *   xi1 (2) - auxiliary field holding the non-time-derivative RHS of the phi
   *             equation (PF3 "Auxiliary")
   *
   * Solve order each increment reproduces PF3: u and phi are updated explicitly
   * from the previous increment's (u, phi, xi1); then xi1 is recomputed from the
   * updated u and phi. At increment 0, u and phi get their initial conditions and
   * xi1 is computed from them.
   */
  std::vector<FieldAttributes> fields = {FieldAttributes("u"),
                                         FieldAttributes("phi"),
                                         FieldAttributes("xi1")};

  SolveBlock explicits;
  explicits.id               = 0;
  explicits.solve_type       = Explicit;
  explicits.solve_timing     = Initialized;
  explicits.field_indices    = {0, 1};
  explicits.dependencies_rhs = make_dependency_set(fields,
                                                   {"old_1(u)",
                                                    "grad(old_1(u))",
                                                    "old_1(phi)",
                                                    "grad(old_1(phi))",
                                                    "old_1(xi1)"});

  SolveBlock xi1_solve;
  xi1_solve.id               = 1;
  xi1_solve.solve_type       = Explicit;
  xi1_solve.solve_timing     = Uninitialized;
  xi1_solve.field_indices    = {2};
  xi1_solve.dependencies_rhs = make_dependency_set(fields, {"u", "phi", "grad(phi)"});

  std::vector<SolveBlock> solves({explicits, xi1_solve});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solves,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
