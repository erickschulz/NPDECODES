- 6-6.f) rename **BoundaryWave.cc** to **boundarywave.cc**
- 6-6.g) rename **BoundaryWave.cc** to **boundarywave.cc**
- 6-6.h) rename **BoundaryWave.h** to **boundarywave.h**
		 replace all occurences of **coord_t** with **Eigen::Vector2d** (probably done automatically)
- 6-6.i) replace in the problem description:
		 **Eigen::VectorXd solveBoundaryWave(
    	 const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fes,
    	 FUNCTOR_U &&u0, FUNCTOR_V &&v0, double T, unsigned int M)**
		 with
		 **Eigen::VectorXd solveBoundaryWave(
         const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p, 
		 FUNCTOR_U &&u0, FUNCTOR_V &&v0, double T, unsigned int N)**

		 
