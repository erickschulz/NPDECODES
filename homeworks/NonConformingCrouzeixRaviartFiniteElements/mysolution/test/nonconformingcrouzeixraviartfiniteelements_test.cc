/**
 * @file
 * @brief NPDE homework NonConformingCrouzeixRaviartFiniteElements code
 * @author Anian Ruoss
 * @date   18.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>

#include <Eigen/Core>

#include "../crfespace.h"
#include "../crl2errordirichletbvp.h"

using namespace NonConformingCrouzeixRaviartFiniteElements;

TEST(CRReferenceFiniteElement, RefEl) {
  CRReferenceFiniteElement cr_ref_el;
  EXPECT_TRUE(cr_ref_el.RefEl() == lf::base::RefEl::kTria());
}

TEST(CRReferenceFiniteElement, Degree) {
  CRReferenceFiniteElement cr_ref_el;
  EXPECT_EQ(cr_ref_el.Degree(), 1);
}

TEST(CRReferenceFiniteElement, NumRefShapeFunctions) {
  CRReferenceFiniteElement cr_ref_el;
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(), 3);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(0), 0);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(1), 1);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(2), 0);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(1, 0), 1);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(1, 1), 1);
  EXPECT_EQ(cr_ref_el.NumRefShapeFunctions(1, 2), 1);
}

TEST(CRReferenceFiniteElement, EvalReferenceShapeFunctions) {
  CRReferenceFiniteElement cr_ref_el;

  Eigen::MatrixXd ref_coords(2, 3);
  ref_coords << 0, 1, 2, 3, 4, 5;

  Eigen::MatrixXd ref_fun_evals(3, 3);
  ref_fun_evals << -5, -7, -9, 5, 9, 13, 1, -1, -3;

  EXPECT_EQ(cr_ref_el.EvalReferenceShapeFunctions(ref_coords), ref_fun_evals);
}

TEST(CRReferenceFiniteElement, GradientsReferenceShapeFunctions) {
  CRReferenceFiniteElement cr_ref_el;

  Eigen::MatrixXd ref_coords(2, 3);
  ref_coords << 0, 1, 2, 3, 4, 5;

  Eigen::MatrixXd ref_fun_grads(3, 6);
  ref_fun_grads << 0, -2, 0, -2, 0, -2, 2, 2, 2, 2, 2, 2, -2, 0, -2, 0, -2, 0;

  EXPECT_EQ(cr_ref_el.GradientsReferenceShapeFunctions(ref_coords),
            ref_fun_grads);
}

TEST(CRReferenceFiniteElement, EvaluationNodes) {
  CRReferenceFiniteElement cr_ref_el;

  Eigen::MatrixXd eval_nodes(2, 3);
  eval_nodes << .5, .5, 0, 0, .5, .5;

  EXPECT_EQ(cr_ref_el.EvaluationNodes(), eval_nodes);
}

TEST(CRReferenceFiniteElement, NumEvaluationNodes) {
  CRReferenceFiniteElement cr_ref_el;
  EXPECT_EQ(cr_ref_el.NumEvaluationNodes(), 3);
}

TEST(CRReferenceFiniteElement, NodalValuesToDofs) {
  CRReferenceFiniteElement cr_ref_el;

  Eigen::MatrixXd nodvals(1, 3);
  nodvals << 0, 1, 2;

  EXPECT_EQ(cr_ref_el.NodalValuesToDofs(nodvals), nodvals);
}

TEST(CRFeSpace, Constructor) {
  std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  CRFeSpace fe_space(mesh_ptr);

  EXPECT_EQ(fe_space.ShapeFunctionLayout(lf::base::RefEl::kSegment()), nullptr);
  EXPECT_EQ(typeid(*(fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()))),
            typeid(CRReferenceFiniteElement));
  EXPECT_EQ(fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()), nullptr);
}

TEST(NonConformingCrouzeixRaviartFiniteElements,
     L2errorCRDiscretizationDirichletBVP) {
  std::vector<double> l2_errors = {0.0227969, 0.00579489, 0.0014546923,
                                   0.000364046};

  // Loop over meshes
  for (int i = 1; i <= 4; ++i) {
    std::string mesh_file = CURRENT_SOURCE_DIR "/../../meshes/refined_square" +
                            std::to_string(i) + ".msh";

    EXPECT_FLOAT_EQ(L2errorCRDiscretizationDirichletBVP(mesh_file),
                    l2_errors[i - 1]);
  }
}
