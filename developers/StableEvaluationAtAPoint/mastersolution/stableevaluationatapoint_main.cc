/**
 * @ file stableevaluationatapoint_main.cc
 * @ brief NPDE homework StableEvaluationAtAPoint
 * @ author Amélie Loher
 * @ date 22.04.20
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>

#include <Eigen/Core>
#include <iostream>
#include <string>

#include "stableevaluationatapoint.h"

using namespace StableEvaluationAtAPoint;

int main(int /*argc*/, const char ** /*argv*/) {
  /* LOADING COARSE MESH */
  // Load mesh into a Lehrfem++ object
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR "/../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();
  // Finite Element Space
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Initial dofhandler
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  lf::base::size_type N_dofs = dofh.NumDofs();

  /* EXACT SOLUTION AND CHOSEN POINT x INSIDE THE DOMAIN*/
  auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };
  Eigen::Vector2d x(0.3, 0.4);

  /* INITIALIZING ERROR ANALYSIS TOOLS AND OBJECTS */
  int N_meshes = 8;  // total number of meshes (coarse + refinement)

  Eigen::VectorXd mesh_sizes(N_meshes);
  mesh_sizes.setZero();
  mesh_sizes(0) = getMeshSize(mesh_p);

  Eigen::VectorXd dofs(N_meshes);
  dofs.setZero();
  dofs(0) = N_dofs;

  // Naive point evaluation
  Eigen::VectorXd errors_Eval(N_meshes);
  errors_Eval.setZero();
  errors_Eval(0) = pointEval(mesh_p);

  // Stable point evaluation
  Eigen::VectorXd errors_stabEval(N_meshes);
  errors_stabEval.setZero();
  Eigen::VectorXd ux(N_meshes);
  ux.setZero();

  /* CONVERGENCE ANALYSIS */
  ux(0) = stab_pointEval(fe_space, u, x);
  errors_stabEval(0) = std::abs(u(x) - ux(0));

  for (int k = 1; k < N_meshes; k++) {  // for each mesh refinement
    // Load finer mesh
    std::string idx = std::to_string(k);
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR
                                                           "/../meshes/square" +
                                                           idx + ".msh");
    mesh_p = reader.mesh();
    fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    lf::base::size_type N_dofs = dofh.NumDofs();

    // Update objects info and evaluate error
    mesh_sizes(k) = getMeshSize(mesh_p);
    dofs(k) = N_dofs;
    errors_Eval(k) = pointEval(mesh_p);

    ux(k) = stab_pointEval(fe_space, u, x);
    errors_stabEval(k) = std::abs(u(x) - ux(k));
    std::cout << "u(x)=" << u(x) << ", ux=" << ux(k) << std::endl;
  }

  // Computing rates of convergence
  double ratesEval[N_meshes - 1];
  double ratesStabEval[N_meshes - 1];
  double log_denum;
  for (int k = 0; k < N_meshes - 1; k++) {
    log_denum = log(mesh_sizes[k] / mesh_sizes[k + 1]);
    ratesEval[k] = log(errors_Eval[k] / errors_Eval[k + 1]) / log_denum;
    ratesStabEval[k] =
        log(errors_stabEval[k] / errors_stabEval[k + 1]) / log_denum;
  }

  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "       Errors for StableEvaluationAtAPoint		"
            << std::endl;
  std::cout << "--------------------- ERRORS ---------------------------"
            << std::endl;
  std::cout << "mesh size"
            << "\t| pointEval"
            << "\t\t| stab_pointEval" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int i = 0; i < 5; i++) {
    std::cout << mesh_sizes(i) << "\t"
              << "\t" << errors_Eval(i) << "\t \t" << errors_stabEval(i)
              << std::endl;
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "      Convergence rates for NAIVE point evaluation       "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| errors_Eval"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << errors_Eval[k];
    if (k > 0) {
      std::cout << "\t\t|" << ratesEval[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  std::cout << "\n" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "      Convergence rates for STABLE point evaluation       "
            << std::endl;
  std::cout << "--------------------- RESULTS ---------------------------"
            << std::endl;
  std::cout << "Iteration"
            << "\t| errors_stabEval"
            << "\t\t| rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k << "\t"
              << "\t|" << errors_stabEval[k];
    if (k > 0) {
      std::cout << "\t\t|" << ratesStabEval[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  // Define output file format
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  std::ofstream file;
  file.open("errors_Eval.csv");
  file << errors_Eval.format(CSVFormat);
  file.close();
  std::cout << "Generated " CURRENT_BINARY_DIR "/errors_Eval.csv" << std::endl;

  return 0;
}
