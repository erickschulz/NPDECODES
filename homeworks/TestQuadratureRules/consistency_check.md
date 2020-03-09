Problem 2-12: Testing built-in quadrature rules of LEHRFEM++

## Remarks R.H.:

Use `EXPECT_` instead of `ASSERT_` in unit tests, because it is desirable that several tests can be conducted instead of aborting execution prematurely.

## Consistency with problem statement

- In 2-12 b): Instead of **boost::math::factorial<double>(I)** we use in the mastersolution **factorial(I)**, where we define the function **factorial()** outside **testQuadOrderTria()**. **(fixed)**

- In 2-12 b) and c): No typedefs in mastersolution (Replace **Vec** and **Mat** with **Eigen::VectorXd** and **Eigen::MatrixXd**, respectively.) **(fixed)**

- In 2-12 b), c) and d): Replace **BOOST_ASSERT** with **assert** (works with #iclude <cassert>) **(fixed)**

- in 2-12 b), c) and d) (PDF): add **in file testquadrules.cc**, and replace qr by quad_rule in function argument name. (const lf::quad::QuadRule &~~qr~~ `quad_rule`)). **(fixed)**

- in 2-12 c): change

  ```c++
  for(int I = 1; I < order; I++){
    for(int J = 0; J < order - I; J++){
  ```

  to

  ```c++
  for(int I = 0; I < order; I++){
    for(int J = 0; J < order; J++){
  ```

  in function `testQuadOrderQuad` in file **testquadrules.cc**. And change corresponding pare in solution.  **(fixed)**

- in **test/testquadrulestest.cc**: remove the printout that prints the whole matrix in `TEST(TestQuadratureRules, TestQuadratureTria)` and `TEST(TestQuadratureRules, TestQuadratureQuad)`. And add another unit test `TEST(TestQuadratureRules, calQuadOrder)` for the function `calcQuadOrder`. **(unit test doesn't show up in PDF anyway)**

