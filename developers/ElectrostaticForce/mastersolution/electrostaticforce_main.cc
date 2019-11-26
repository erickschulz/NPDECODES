/**
 * @file newproblem_main.cc
 * @brief NPDE homework NewProblem code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

#include "electrostaticforce.h"

using namespace ElectrostaticForce;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), CURRENT_SOURCE_DIR
                                  "/../meshes/emforce3.msh");
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>

  // Finite element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Solve BVP
  Eigen::VectorXd discrete_solution = solvePoissonBVP(fe_space_p);


  /* Output results to vtk file */
  // We store data by keeping only the coefficients of nodal basis functions
  // In that sense, we are plotting the values of the solution at the vertices
  lf::io::VtkWriter vtk_writer(mesh_p, "ElectrostaticForcePoissonBVP_solution.vtk");
  // Write nodal data taking the values of the discrete solution at the vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    if (dofh.Entity(global_idx).RefEl() == lf::base::RefElType::kPoint) {
      nodal_data->operator()(dofh.Entity(global_idx)) = discrete_solution[global_idx];
    }
  };
  vtk_writer.WritePointData("ElectrostaticForcePoissonBVP_solution", *nodal_data);
  /* SAM_LISTING_END_1 */
  std::cout << "\n The solution vector was written to:" << std::endl;
  std::cout << ">> ElectrostaticForcePoissonBVP_solution.vtk\n" << std::endl;

  return 0;
}
