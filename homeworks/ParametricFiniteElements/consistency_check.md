- Replace in all problem statements **size_t** with **unsigned int**
  (occurs in (5.6e), (5.6f), (5.6h),(5.6i), (5.6j), (5.6k), (5.6l) ) 

- Replace in all problem statements: the file for the solution is now **parametricfiniteelements.h** instead of **main.cpp** !!

- (5.6e): Replace **template <typename Functor>** with **template <typename FUNCTOR>** and correspondingly **Functor Psi** with **FUNCTOR &&Psi**

- (5.6f): Replace **template <typename Functor1, typename Functor2>** with **template <typename FUNCTOR1, typename FUNCTOR2>** and correspondingly **Functor1 Psi, Functor2 alpha** with **FUNCTOR1 &&Psi, FUNCTOR2 &&alpha**

- (5.6f) Hint 1: the shape functions and their gradients on the reference elements are provided by the functions
		 
		 Eigen::Vector4d bhats(Eigen::Vector2d xhat);

		 and 

		 Eigen::MatrixXd bhats_grad(Eigen::Vector2d xhat);

- (5.6h): replace **int geoThermLocal2Global()** with **int geoThermLocalToGlobal()**

- (5.6i): Replace **template <typename Functor1, typename Functor2>** with **template <typename FUNCTOR1, typename FUNCTOR2>** and correspondingly **Functor1 Psi, Functor2 alpha** with **FUNCTOR1 &&Psi, FUNCTOR2 &&alpha**

		  Replace **std::vector<Triplets> assemGeoThermMat()** with

		  **std::vector<Eigen::Triplet<double>> assembleGeoTherm()**

- (5.6j): Replace **std::vector<Triplets>** with **std::vector<Eigen::Triplet<double>>**

		 Remove **template <typename Functor1, typename Functor2>**


- (5.6k): Replace **template <typename Functor1, typename Functor2>** with **template <typename FUNCTOR1, typename FUNCTOR2>** and correspondingly **Functor1 Psi, Functor2 alpha** with **FUNCTOR1 &&Psi, FUNCTOR2 &&alpha**



- (5.6l): Replace **template <typename Functor>** with **template <typename FUNCTOR>** and correspondingly **Functor Psi** with **FUNCTOR &&Psi**




I assume that the codes will be updated anyways.
