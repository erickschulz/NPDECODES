/**
 * @ file stableevaluationatapoint_main.cc 
 * @ brief NPDE homework StableEvaluationAtAPoint 
 * @ author Amélie Loher
 * @ date 22.04.20
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "stableevaluationatapoint.h"

#include <iostream>
#include <string>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>

using namespace StableEvaluationAtAPoint;

int main(int /*argc*/, const char** /*argv*/) {
  
  // Load mesh into a Lehrfem++ object
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init), CURRENT_SOURCE_DIR "/../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();
  // Finite Element Space
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space 
  									= std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Initial dofhandler
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  lf::base::size_type N_dofs = dofh.NumDofs();

  Eigen::VectorXd mesh_sizes(5); mesh_sizes.setZero();
  Eigen::VectorXd dofs(5); dofs.setZero();
  Eigen::VectorXd errors_Eval(5); errors_Eval.setZero();
  Eigen::VectorXd errors_stabEval(5); errors_stabEval.setZero();

  Eigen::VectorXd ux(5); ux.setZero();

  mesh_sizes(0) = getMeshSize(mesh_p);
  dofs(0) = N_dofs;
  errors_Eval(0) = pointEval(mesh_p);

  auto u = [] (Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log( (x + one).norm() );
  };
  Eigen::Vector2d x(0.3, 0.4);
  
  ux(0) = stab_pointEval(fe_space, u, x);
  errors_stabEval(0) = std::abs(u(x) - ux(0));
  
  	
  for(int k = 1; k < 5; k++) {
    std::string idx = std::to_string(k);
    
	auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
	lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR "/../meshes/square" + idx + ".msh");
	mesh_p = reader.mesh();
    fe_space = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    lf::base::size_type N_dofs = dofh.NumDofs();

    mesh_sizes(k) = getMeshSize(mesh_p);
    dofs(k) = N_dofs;
    errors_Eval(k) = pointEval(mesh_p);

    ux(k) = stab_pointEval(fe_space, u, x);
    errors_stabEval(k) = std::abs(u(x) - ux(k));
    std::cout << "u(x)" << u(x) << ", ux" << ux(k) << std::endl;
  }


  std::cout << "*********************************************************" << std::endl;
  std::cout << "       ERRORS FOR STABLEEVALUATIONATAPOINT				 " << std::endl;
  std::cout << "*********************************************************" << std::endl;
  
  std::cout << "---------------------------------------------------------" << std::endl;
  std::cout << "mesh size"
  			<< "\t| error pointEval"
			<< "\t\t| error stab_pointEval" << std::endl;


  for(int i = 0; i < 5; i++) {
    std::cout << mesh_sizes(i) << "\t"
			  << "          " << errors_Eval(i) << "\t \t"
    		  << "          " << errors_stabEval(i) << std::endl;
  }
  std::cout << "---------------------------------------------------------" << std::endl;
 
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
