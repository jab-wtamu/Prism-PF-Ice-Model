// =====================================================================================
// customPDE.h  (ICE APP: u, phi)
// =====================================================================================
//
// Changes in this version
//
// - Added model constant:
//     Dtilde  (reduced diffusion coefficient for u equation)
//     lambda (for phi source coupling term)

// Everything else remains minimal to keep the app stable.
// =====================================================================================

#include <core/matrixFreePDE.h>

using namespace dealii;

template <int dim, int degree>
class customPDE : public MatrixFreePDE<dim, degree>
{
public:
  customPDE(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
    , userInputs(_userInputs)
  {}

  void
  setInitialCondition([[maybe_unused]] const Point<dim>  &p,
                      [[maybe_unused]] const unsigned int index,
                      [[maybe_unused]] double            &scalar_IC,
                      [[maybe_unused]] Vector<double>    &vector_IC) override;
//Declares: you override the base method that sets initial conditions.

  void
  setNonUniformDirichletBCs([[maybe_unused]] const Point<dim>  &p,
                            [[maybe_unused]] const unsigned int index,
                            [[maybe_unused]] const unsigned int direction,
                            [[maybe_unused]] const double       time,
                            [[maybe_unused]] double            &scalar_BC,
                            [[maybe_unused]] Vector<double>    &vector_BC) override;
//Declares: you override boundary condition assignment.
//original*
private:
#include <core/typeDefs.h>

  const userInputParameters<dim> userInputs;

  void
  explicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  void
  nonExplicitEquationRHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

  void
  equationLHS(
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;

// Function to set postprocessing expressions (in postprocess.h)
#ifdef POSTPROCESS_FILE_EXISTS
  void
  postProcessedFields(
    [[maybe_unused]] const variableContainer<dim, degree, VectorizedArray<double>> &variable_list,
    [[maybe_unused]] variableContainer<dim, degree, VectorizedArray<double>>       &pp_variable_list,
    [[maybe_unused]] const Point<dim, VectorizedArray<double>>                     q_point_loc,
    [[maybe_unused]] const VectorizedArray<double> element_volume) const override;
#endif
//original*

#ifdef NUCLEATION_FILE_EXISTS
  double
  getNucleationProbability([[maybe_unused]] variableValueContainer variable_value,
                           [[maybe_unused]] double                 dV) const override;
#endif
  //removed model constants for now
  // ================================================================
  // Model constants needed by IC scaffold (phi seeds)
  // ================================================================
  Tensor<1, dim> center1 = userInputs.get_model_constant_rank_1_tensor("center1");
  Tensor<1, dim> center2 = userInputs.get_model_constant_rank_1_tensor("center2");
  //Loads two vector-like parameters (center1, center2) from the input file.
  //Used to place seed centers.
  double radius1 = userInputs.get_model_constant_double("radius1");
  double radius2 = userInputs.get_model_constant_double("radius2");
  //Loads radii of the seeds.
  // ================================================================
  // NEW: diffusion coefficient for u equation
  // ================================================================
  double Dtilde = userInputs.get_model_constant_double("Dtilde");
  
  // Coupling constant for phi equation 
  double lambda = userInputs.get_model_constant_double("lambda");

};
