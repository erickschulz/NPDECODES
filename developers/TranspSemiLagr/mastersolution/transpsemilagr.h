#include<lf/uscalfe/uscalfe.h>
#include <Eigen/Core>
#include <Eigen/Sparse>


namespace TranspSemiLagr{


void enforce_zero_boundary_conditions(std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space, lf::assemble::COOMatrix<double> &A, Eigen::VectorXd &b){
  
  lf::mesh::utils::MeshFunctionGlobal mf_zero{[](const Eigen::Vector2d& /*x*/){return 0.0;}};
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto flag_values{lf::uscalfe::InitEssentialConditionFromFunction(
        fe_space->LocGlobMap(), *rsf_edge_p, bd_flags,mf_zero)};

   lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&flag_values](lf::assemble::glb_idx_t dof_idx) {
        return flag_values[dof_idx];
      },
      A, b);
}



template<typename FUNCTOR>
Eigen::VectorXd semiLagr_step( std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
                                          const Eigen::VectorXd & u0_vector,
                                          FUNCTOR v, double tau){

  //Assemble left hand side A = A_lm + tau*A_s
  lf::assemble::COOMatrix<double> A (fe_space->LocGlobMap().NumDofs(), fe_space->LocGlobMap().NumDofs());

  //stiffness matrix 
  lf::uscalfe::ReactionDiffusionElementMatrixProvider stiffness_element_matrix_provider(fe_space, lf::mesh::utils::MeshFunctionConstant(tau),
                                                                      lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(), stiffness_element_matrix_provider, A);

  //lumped mass matrix
  MassLumpedElementMatrixProvider  lumped_mass_element_matrix_provider ([](Eigen::Vector2d /*x*/){return 1.0;});
  lf::assemble::AssembleMatrixLocally(0,fe_space->LocGlobMap(), fe_space->LocGlobMap(),stiffness_element_matrix_provider, A);


  //warp u0 into a mesh function & assemble rhs.
  auto u0_mf = lf::uscalfe::MeshFunctionFE(fe_space, u0_vector);
  UpwindLagrangianElementVectorProvider vector_provider (v, tau, fe_space->Mesh(), u0_mf);
  Eigen::VectorXd b (fe_space->LocGlobMap().NumDofs());
  b.setZero();
  lf::assemble::AssembleVectorLocally(0, fe_space->LocGlobMap(), vector_provider, b);

  
  enforce_zero_boundary_conditions(fe_space,A,b);
  
  //solve LSE
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  return solver.solve(b);
  
}


template<typename FUNCTOR>
Eigen::VectorXd solverot( 
         std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
         FUNCTOR u0,
         int N, double T){

  double tau = T / N;

  Eigen::VectorXd u0_vector  = lf::uscalfe::NodalProjection(*fe_space, lf::mesh::utils::MeshFunctionGlobal(u0));
  auto v = [](Eigen::Vector2d x){return (Eigen::Vector2d()  << -x(1) + 3.0*x(0)*x(0), x(0)).finished();};

  for(int i = 0; i < N; ++i){
    u0_vector = semiLagr_step(fe_space, u0_vector, v, tau);
  }

  return u0_vector;
}

template<typename FUNCTOR>
Eigen::VectorXd reaction_step(std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space, 
                                          const Eigen::VectorXd& u0_vector,FUNCTOR c, double tau){
   MassLumpedElementMatrixProvider  mass_element_matrix_provider_c (c);
   MassLumpedElementMatrixProvider  mass_element_matrix_provider_1 ([](Eigen::Vector2d /*x*/){return 1.0;});

   lf::assemble::COOMatrix<double> mass_matrix_c (fe_space->LocGlobMap().NumDofs(), fe_space->LocGlobMap().NumDofs());
   lf::assemble::COOMatrix<double> mass_matrix_1 (fe_space->LocGlobMap().NumDofs(), fe_space->LocGlobMap().NumDofs());

   lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(), mass_element_matrix_provider_c,mass_matrix_c);
   lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(),  mass_element_matrix_provider_1,mass_matrix_1);

   Eigen::SparseMatrix<double> mass_matrix_c_sparse = mass_matrix_c.makeSparse();
   Eigen::SparseMatrix<double> mass_matrix_1_sparse = mass_matrix_1.makeSparse();
   
   Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
   solver.compute(mass_matrix_1_sparse);

   Eigen::VectorXd k1 = solver.solve(mass_matrix_c_sparse * u0_vector);
   Eigen::VectorXd k2 = solver.solve(mass_matrix_c_sparse* (u0_vector + 0.5*tau*k1));

  return u0_vector + tau*k2;
 }


template<typename FUNCTOR>
Eigen::VectorXd solvetrp(std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
                                             FUNCTOR u0, int N, double T){
  double tau = T/N;

  Eigen::VectorXd u0_vector =  lf::uscalfe::NodalProjection(*fe_space, lf::mesh::utils::MeshFunctionGlobal(u0));

  auto v = [](Eigen::Vector2d x){return (Eigen::Vector2d() << -x(1) + 3.0*x(0)*x(0), x(0)).finished();};
  auto c = [](Eigen::Vector2d x){return -6.0*x(0);};

  //Strang splitting scheme
  //first SemiLagr half step:
  u0_vector = semiLagr_step(fe_space, u0_vector, v, 0.5*tau);

  //intermediate time steps: Combine two semiLagr half steps to one step
  for(int i = 0; i < N-1; ++i){
     u0_vector = reaction_step(fe_space, u0_vector, c, tau);
     u0_vector = semiLagr_step(fe_space, u0_vector, v, tau);
  }

  //final reaction step and semiLagr half step
  u0_vector = reaction_step(fe_space, u0_vector, c, tau);
  u0_vector = semiLagr_step(fe_space,u0_vector, v, tau*0.5);

  return u0_vector;
}

} //namespace TranspSemiLagr
