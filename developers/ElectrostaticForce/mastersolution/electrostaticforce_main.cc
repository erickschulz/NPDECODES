/**
 * @file electrostaticforce_main.cc
 * @brief ElectrostaticForce
 * @author Erick Schulz
 * @date 27.11.2019
 * @copyright Developed at ETH Zurich
 */

#include "electrostaticforce.h"

using namespace ElectrostaticForce;

// Spcify output format for table
const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  std::cout << "\nPROBLEM - ElectrostaticForce\n" << std::endl;

  // TOOLS AND DATA
  int N_meshes = 5;            // No. of meshes for convergence test
  Eigen::VectorXd approx_sol;  // basis coeff expansion of approx solution
  Eigen::Vector2d approx_force_boundary_functional;
  Eigen::Vector2d approx_force_domain_functional;
  double errorsL2PoissonBVP[N_meshes];  // L2 error domain BVP
  double
      errorsl2ForceBoundaryFunctional[N_meshes];   // l2 err force bd functional
  double errorsl2ForceDomainFunctional[N_meshes];  // l2 err force domain
                                                   // functional
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p;
  std::shared_ptr<const lf::mesh::Mesh> mesh_p;
  double mesh_sizes[N_meshes];  // meshwidths

  // CREATE "EXACT" SOLUTION AND FORCE VECTOR
  // Create analytic (exact) solution of the Poisson BVP as mesh function
  // See sub-problem b)
  /* SAM_LISTING_BEGIN_A */
  Eigen::Vector2d a(-16.0 / 15.0, 0.0);
  Eigen::Vector2d b(-1.0 / 15.0, 0.0);
  // Lambda function providing exact solution
  auto uExact = [&a, &b](Eigen::VectorXd x) -> double {
    return (log((x - a).norm()) - log((x - b).norm())) / log(2) - 1;
  };
  /* SAM_LISTING_END_A */
  // Mesh function wrapper for exact solution
  lf::mesh::utils::MeshFunctionGlobal mf_uExact{uExact};
  // Compute "exact" force using overkill quadrature
  Eigen::Vector2d exact_force = computeExactForce();

  for (int i = 0; i < N_meshes; i++) {  // for each mesh
    // READ MESH INTO LEHRFEMPP
    // Load mesh into a Lehrfem++ object
    std::string mesh_file = CURRENT_SOURCE_DIR "/../meshes/emforce" +
                            std::to_string(i + 1) + ".msh";
    auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
    mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>
    mesh_sizes[i] = getMeshSize(mesh_p);
    // Finite element space
    fe_space_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

    // SOLVE POISSON DIRICHLET BVP
    approx_sol = solvePoissonBVP(fe_space_p);
    // COMPUTE FORCE FUNCTIONALs
    approx_force_boundary_functional =
        computeForceBoundaryFunctional(fe_space_p, approx_sol);
    approx_force_domain_functional =
        computeForceDomainFunctional(fe_space_p, approx_sol);

    /* SAM_LISTING_BEGIN_1 */
    // COMPUTE L2 ERROR FOR THE POISSON DIRICHLET PROBLEM
    auto mf_approx_sol = lf::uscalfe::MeshFunctionFE(fe_space_p, approx_sol);
    errorsL2PoissonBVP[i] = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
        *mesh_p, lf::uscalfe::squaredNorm(mf_uExact - mf_approx_sol), 2));
    /* SAM_LISTING_END_1 */

    // COMPUTE l2 ERROR FOR THE FORCE
    errorsl2ForceBoundaryFunctional[i] =
        (exact_force - approx_force_boundary_functional).norm();
    errorsl2ForceDomainFunctional[i] =
        (exact_force - approx_force_domain_functional).norm();
  }

  // CALCULATING CONVERGENCE RATES
  double ratesL2PoissonBVP[N_meshes - 1];
  double ratesl2ForceBoundaryFunctional[N_meshes - 1];
  double ratesl2ForceDomainFunctional[N_meshes - 1];
  double log_mesh_size_denumerator;
  for (int k = 0; k < N_meshes - 1; k++) {
    ratesL2PoissonBVP[k] =
        log(errorsL2PoissonBVP[k] / errorsL2PoissonBVP[k + 1]) /
        log(mesh_sizes[k] / mesh_sizes[k + 1]);
    ratesl2ForceBoundaryFunctional[k] =
        log(errorsl2ForceBoundaryFunctional[k] /
            errorsl2ForceBoundaryFunctional[k + 1]) /
        log(mesh_sizes[k] / mesh_sizes[k + 1]);
    ratesl2ForceDomainFunctional[k] =
        log(errorsl2ForceDomainFunctional[k] /
            errorsl2ForceDomainFunctional[k + 1]) /
        log(mesh_sizes[k] / mesh_sizes[k + 1]);
  }

  // DISPLAY CONVERGENCE TABLE FOR POISSON PROBLEM
  std::cout.precision(5);
  std::cout << "\n" << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "      ERRORS AND CONVERGENCE RATES FOR PROBLEM        "
            << std::endl;
  std::cout << "               ELECTROSTATIC FORCE                    "
            << std::endl;
  std::cout << "*********************************************************"
            << std::endl;
  std::cout << "" << std::endl;
  std::cout << "    L2 errors of approximate scalar solutions to the     "
            << std::endl;
  std::cout << "   Poisson Problem with essential boundary conditions    "
            << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "iteration"
            << "\t mesh size"
            << "\t L2 error"
            << "\t rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k + 1 << "\t\t|" << std::fixed << mesh_sizes[k]
              << std::scientific << "\t|" << errorsL2PoissonBVP[k];
    if (k > 0) {
      std::cout << std::fixed << "\t|" << ratesL2PoissonBVP[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "" << std::endl;

  // DISPLAY CONVERGENCE TABLE FOR BOUNDARY FORCE FUNCTIONAL
  std::cout << "         l2 errors of the approximate forces             "
            << std::endl;
  std::cout << "          computed with BOUNDARY functional              "
            << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "iteration"
            << "\t mesh size"
            << "\t l2 error"
            << "\t rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k + 1 << "\t\t|" << std::fixed
              << mesh_sizes[k]  //<< std::scientific
              << "\t|" << errorsl2ForceBoundaryFunctional[k];
    if (k > 0) {
      std::cout << std::fixed << "\t|" << ratesl2ForceBoundaryFunctional[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl
            << std::endl;

  // DISPLAY CONVERGENCE TABLE FOR DOMAIN FORCE FUNCTIONAL
  std::cout << "         l2 errors of the approximate forces             "
            << std::endl;
  std::cout << "          computed with DOMAIN functional                "
            << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  std::cout << "iteration"
            << "\t mesh size"
            << "\t l2 error"
            << "\t rates" << std::endl;
  std::cout << "---------------------------------------------------------"
            << std::endl;
  for (int k = 0; k < N_meshes; k++) {
    std::cout << k + 1 << "\t\t|" << std::fixed
              << mesh_sizes[k]  //<< std::scientific
              << "\t|" << errorsl2ForceDomainFunctional[k];
    if (k > 0) {
      std::cout << std::fixed << "\t|" << ratesl2ForceDomainFunctional[k - 1];
    }
    std::cout << "\n";
  }
  std::cout << "---------------------------------------------------------"
            << std::endl;

  /* Output results to vtk file */
  // We store data by keeping only the coefficients of nodal basis functions
  // In that sense, we are plotting the values of the solution at the vertices
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  lf::io::VtkWriter vtk_writer(
      mesh_p, CURRENT_BINARY_DIR "/ElectrostaticForcePoissonBVP_solution.vtk");
  // Write nodal data taking the values of the discrete solution at vertices
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
  for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
    if (dofh.Entity(global_idx).RefEl() == lf::base::RefElType::kPoint) {
      nodal_data->operator()(dofh.Entity(global_idx)) = approx_sol(global_idx);
    }
  };
  vtk_writer.WritePointData("ElectrostaticForcePoissonBVP_solution",
                            *nodal_data);

  std::cout << "\nThe solution vector was written to:" << std::endl;
  std::cout << "ElectrostaticForcePoissonBVP_solution.vtk\n" << std::endl;

  return 0;
}
