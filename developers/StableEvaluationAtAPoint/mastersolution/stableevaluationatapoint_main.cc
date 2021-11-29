/**
 * @file stableevaluationatapoint_main.cc
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher & Erick Schulz
 * @date 22.04.20
 * @copyright Developed at SAM, ETH Zurich
 */

#include "stableevaluationatapoint.h"

#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>

using namespace StableEvaluationAtAPoint;

int main(int /*argc*/, const char ** /*argv*/) {

  //exact solution
  auto uExact = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };
  // Fixed evaluation point
  Eigen::Vector2d x(0.3, 0.4);
  std::cout << "Exact evaluation at (0.3,0.4) : " << uExact(x) << std::endl;

  //Number of meshes used in the error analysis:
  int N_meshes = 7;

  /* INITIALIZING ERROR ANALYSIS TOOLS AND OBJECTS */
  Eigen::VectorXd mesh_sizes(N_meshes);
  mesh_sizes.setZero();

  Eigen::VectorXd dofs(N_meshes);
  dofs.setZero();

  //Error vector used for the error analysis in exercise b)
  Eigen::VectorXd errors_potential(N_meshes);
  errors_potential.setZero();

  //Error vectors used for the error analysis in exercise g/h
  Eigen::VectorXd errors_direct_eval(N_meshes);
  errors_direct_eval.setZero();
  Eigen::VectorXd errors_stable_eval(N_meshes);
  errors_stable_eval.setZero();


  //iterate over meshes:
  for (int k = 0; k < N_meshes; k++) { 
    //read mesh::
    std::string idx = std::to_string(k);
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR
                                                           "/../meshes/square" +
                                                           idx + ".msh");
    auto mesh_p = reader.mesh();
  
    //Initialize Fe-Space and dofh
    auto fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    dofs(k) = dofh.NumDofs();

    // Printing mesh statistics
    mesh_sizes(k) = MeshSize(mesh_p);
    std::cout << "square.msh: "
              << "N_dofs = " << dofs(k) << ", h=" << mesh_sizes(k) << std::endl;

    //error anlysis part b) (Potentials)
    errors_potential(k) = PointEval(mesh_p);

    //error analysis part g/h: Compare direct vs stable point evaluation:
    auto [direct_eval,stable_eval] = ComparePointEval(fe_space,uExact,x);
    errors_direct_eval(k) = std::abs(uExact(x) - direct_eval);
    errors_stable_eval(k) = std::abs(uExact(x) - stable_eval);
  }

  //Compute rates of convergence:
  Eigen::VectorXd rates_potential(N_meshes-1);
  Eigen::VectorXd rates_direct(N_meshes-1);
  Eigen::VectorXd rates_stable(N_meshes-1);

  for(int k = 0; k < N_meshes-1; ++k){
    double log_denum = std::log(mesh_sizes(k)/mesh_sizes(k+1));
    rates_potential(k) = std::log(errors_potential(k)/errors_potential(k+1))/log_denum;
    rates_direct(k) = std::log(errors_direct_eval(k)/errors_direct_eval(k+1))/log_denum;
    rates_stable(k) = std::log(errors_stable_eval(k)/errors_stable_eval(k+1))/log_denum;
  }



  std::cout << "Subtask b) Evaluation based on Potentials \n";
  std::cout << "Errors: " << errors_potential << "\n";
  std::cout << "Rates: " << rates_potential << "\n";
  std::cout << "Subtask h) Comparison of direct and stable evaluations: \n";
  std::cout << "Errors direct: " << errors_direct_eval << "\n";
  std::cout << "Rates direct: " << rates_direct << "\n";
  std::cout << "Errors stable: " << errors_stable_eval << "\n";
  std::cout << "Rates stable: " << rates_stable << "\n";


  // Define output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  Eigen::MatrixXd convergence_potential(N_meshes,2);
  convergence_potential << mesh_sizes, errors_potential;

  Eigen::MatrixXd convergence_stable(N_meshes,3);
  convergence_stable << mesh_sizes, errors_direct_eval, errors_stable_eval;

  std::ofstream file;
  file.open("convergence_potential.csv");
  file << "h, Error u(x) (Potential) \n";
  file << convergence_potential.format(CSVFormat);
  file.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/convergence_potential.csv" << std::endl;


  file.open("convergence_stable.csv");
  file << "h, Error u(x) (Direct), Error u(x) (Stable) \n";
  file << convergence_stable.format(CSVFormat);
  file.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/convergence_stable.csv" << std::endl;

  // Plot
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_convergence_potential.py " CURRENT_BINARY_DIR);
  std::system("python3 " CURRENT_SOURCE_DIR
              "/plot_convergence_stable.py " CURRENT_BINARY_DIR);

  return 0;
}
