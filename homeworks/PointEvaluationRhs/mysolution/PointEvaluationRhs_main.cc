/**
 * @ file PointEvaluationRhs_main.cc
 * @ brief NPDE homework PointEvaluationRhs code
 * @ author Christian Mitsch
 * @ date 22.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <mgl2/mgl.h>
#include <iomanip>
#include "norms.h"
#include "pointEvaluation.h"

void plot_result(const std::vector<double> &dof_a,
                 const std::vector<double> &l2_a,
                 const std::vector<double> &h1_a);

int main() {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(1, 1.0);

  // Start of numerical experiment
  std::vector<double> dof_a{};
  std::vector<double> l2_a{};
  std::vector<double> h1_a{};
  Eigen::VectorXd sol_vec;

  // Runs with initial mesh
  lf::assemble::UniformFEDofHandler dofh_initial(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  auto result = PointEvaluationRhs::normsSolutionPointLoadDirichletBVP(
      dofh_initial, Eigen::Vector2d(1.3, 1.7), sol_vec);
  unsigned N_dofs = dofh_initial.NoDofs();
  dof_a.push_back(N_dofs);
  l2_a.push_back(result.first);
  h1_a.push_back(result.second);

  // Necessary for regular refinement
  auto mesh_factory2 = std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::refinement::MeshHierarchy my_hierarchy{mesh_p, mesh_factory2};

  for (int k = 1; k < 7; k++) {
    my_hierarchy.RefineRegular();

    mesh_p = my_hierarchy.getMesh(k);
    lf::assemble::UniformFEDofHandler dofh(mesh_p,
                                           {{lf::base::RefEl::kPoint(), 1}});
    unsigned N_dofs = dofh.NoDofs();
    dof_a.push_back(N_dofs);
    sol_vec.resize(N_dofs);

    result =
        PointEvaluationRhs::normsSolutionPointLoadDirichletBVP(dofh, Eigen::Vector2d(1.3, 1.7),sol_vec);
    l2_a.push_back(result.first);
    h1_a.push_back(result.second);
    // Write vtk file
    std::stringstream filename;
    filename << "rhseval" << k << ".vtk";
    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    // need the newest pointer
    auto mds = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (auto &node : mesh_p->Entities(2)) {
      mds->operator()(node) = sol_vec(dofh.GlobalDofIndices(node)[0]);
    }
    vtk_writer.WritePointData("solution_data", *mds);
  }

  // Plot the resulting norms
  plot_result(dof_a, l2_a, h1_a);

  // Print to std output
  std::cout << " dof      l2         h1 " << std::endl;
  for (int i = 0; i < dof_a.size(); i++) {
    std::cout << std::setw(5) << dof_a.at(i) << "   " << std::setw(5)
              << l2_a.at(i) << "   " << std::setw(5) << h1_a.at(i) << std::endl;
  }

  // END_SOLUTION
}

void plot_result(const std::vector<double> &dof_a,
                 const std::vector<double> &l2_a,
                 const std::vector<double> &h1_a) {
  mglGraph graph;
  mglData x_(dof_a.size(), dof_a.data());
  mglData y_(l2_a.size(), l2_a.data());
  mglData z_(h1_a.size(), h1_a.data());
  graph.SetRange('y', 0, 1);
  graph.SetRange('x', 0, *std::max_element(dof_a.begin(), dof_a.end()));

  graph.Title("Convergence");
  graph.Axis();
  graph.Label('x', "Ndofs", 0);
  graph.Label('y', "Norms", 0);
  graph.Plot(x_, y_);
  graph.Plot(x_, z_);

  graph.WriteFrame("convergence.png");
}
