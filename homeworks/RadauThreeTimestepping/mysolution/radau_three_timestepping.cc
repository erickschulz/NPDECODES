/**
 * @file
 * @brief NPDE homework RadauThreeTimestepping
 * @author Erick Schulz
 * @date 08/04/2019
 * @copyright Developed at ETH Zurich
 */

#include "radau_three_timestepping.h"

namespace RadauThreeTimestepping {

/** @brief Implementation of the right hand side (time dependent) source vector
 * for the parabolic heat equation*/
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd rhsVectorheatSource(const lf::assemble::DofHandler &dofh,
                                    double time) {
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NoDofs());
  // Right-hand side vector has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
    /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return phi;
}
/* SAM_LISTING_END_1 */

/** @Brief Heat evolution solver: the solver obtains the
discrete evolution operator from the Radau3MOLTimestepper class and repeatedly
iterates its applicaiton starting from the initial condition
* @param m is total number of steps until final time final_time (double) */
/* SAM_LISTING_BEGIN_6 */
Eigen::VectorXd solveHeatEvolution(const lf::assemble::DofHandler &dofh,
                                   unsigned int m, double final_time) {
  Eigen::VectorXd discrete_heat_sol;
   /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return discrete_heat_sol;
}
/* SAM_LISTING_END_6 */

/* Implementing member function Eval of class LinFEMassMatrixProvider*/
Eigen::Matrix<double, 3, 3> LinFEMassMatrixProvider::Eval(
    const lf::mesh::Entity &tria) {
  Eigen::Matrix<double, 3, 3> elMat;
    /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return elMat;  // return the local mass element matrix
}  // LinFEMassMatrixProvider::Eval

/* Implementing constructor of class Radau3MOLTimestepper */
/* SAM_LISTING_BEGIN_4 */
Radau3MOLTimestepper::Radau3MOLTimestepper(const lf::assemble::DofHandler &dofh)
    : dofh_(dofh) {
     /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
}  // Radau3MOLTimestepper constructor
/* SAM_LISTING_END_4 */

/* Implementation of Radau3MOLTimestepper member functions */
// The function discreteEvolutionOperator() returns the discretized evolution
// operator as obtained from the Runge-Kutta Radau IIA 2-stages method using the
// Butcher table as stored in the Radau3MOLTimestepper class
/* SAM_LISTING_BEGIN_5 */
Eigen::VectorXd Radau3MOLTimestepper::discreteEvolutionOperator(
    double time, double tau, const Eigen::VectorXd &mu) const {
  Eigen::VectorXd discrete_evolution_operator;
     /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return discrete_evolution_operator;
}
/* SAM_LISTING_END_5 */

}  // namespace RadauThreeTimestepping
