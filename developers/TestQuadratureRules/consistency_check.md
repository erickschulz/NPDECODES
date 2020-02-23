Problem 2-12: Testing built-in quadrature rules of LEHRFEM++

- In 2-12 b): Instead of **boost::math::factorial<double>(I)** we use in the mastersolution **factorial(I)**, where we define the function **factorial()** outside **testQuadOrderTria()**.

- In 2-12 b) and c): No typedefs in mastersolution (Replace **Vec** and **Mat** with **Eigen::VectorXd** and **Eigen::MatrixXd**, respectively.)

- In 2-12 b), c) and d): Replace **BOOST_ASSERT** with **assert** (works with #iclude <cassert>)




