
- (3-11.a): 
	- The two functions   G and gradG are now implemented in a struct
		
			class FundamentalSolution {
			public:
			FundamentalSolution(Eigen::Vector2d x) : x_{x} {}

			// Computes  G_x(y)
			double operator()(Eigen::Vector2d y);
			// Computes the gradient of  G_x(y)
			Eigen::Vector2d grad(Eigen::Vector2d y);

			private:
			Eigen::Vector2d x_;
			}; 

	- Mention that the shortcut to evaluate the normals on the unit square is implemented in the function 
	    
			Eigen::Vector2d OuterNormalUnitSquare(Eigen::Vector2d x);

- (3-11.b): 
	- Function signautre was updated to 

        	double PointEval(std::shared_ptr<const lf::mesh::Mesh> mesh_p);
	
	- New solution plot using python script **plot_convergence_potential.py**



- (3-11.f): 
	- Functions are now also collected in a small struct:

			class Psi {
				public:
				Psi(Eigen::Vector2d center) : center_(center) {}

				// computes Psi_x(y)
				double operator()(Eigen::Vector2d y);
				// computes grad(Psi_x)(y)
				Eigen::Vector2d grad(Eigen::Vector2d y);
				// computes the laplacian of Psi_x at y
				double lapl(Eigen::Vector2d y);

				private:
				Eigen::Vector2d center_ = Eigen::Vector2d(0.5, 0.5);
			};

	- Function signature was updated to

			double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
				Eigen::VectorXd uFE, const Eigen::Vector2d x);

	- Mention that the following function additionally checks the assumptions on the point x and then evaluates Jstar

			double StablePointEvaluation(
				std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
				Eigen::VectorXd uFE, const Eigen::Vector2d x);

	


- (3-11.g): 
	- Function signature of the Laplace solver was updated to
			
			template <typename FUNCTOR>
			Eigen::VectorXd SolveBVP(
				const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
				FUNCTOR &&g)

	- Instead of **Jstar**, the function **StablePointEvaluation** should be used to perform the stable point evalation (checks the assumptions on x)
 
	- For the direct evaluation, mention the following function to evaluate the FE solution at the global coordinates x.

			template <typename FUNCTOR>
			std::pair<double, double> ComparePointEval(
				std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
				FUNCTOR &&g, Eigen::Vector2d x)

- (3-11.h): - New convergence plot based on the python script **plot_convergence_stable.py**. The plot now also shows the convergence for the direct evaluation method.


