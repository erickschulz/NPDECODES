- Remove in the problem description any occurrences of **BETL**
- Adapt all function signatures accordingly

- (3-11.a): Replace the current function signature with

			**template <typename FUNCTOR>
			double PSL(std::shared_ptr<const lf::mesh::Mesh> mesh, const FUNCTOR &v,
           	const Eigen::Vector2d x)**

		   	and

			**template <typename FUNCTOR>
			double PDL(std::shared_ptr<const lf::mesh::Mesh> mesh, const FUNCTOR &v,
           	const Eigen::Vector2d x)**

			Adapt the description of parameters (input arguments).

- (3-11.b): Remove matlab script. Instead, they shall modify the main function **stableevaluationatapoint_main.cc** accordingly. 

- (3-11.f): Remove function signature **numeric_t Psi(const coord_t& y, coord_t& gradPsi, numeric_t& laplPsi);**

			Instead explain that these functions are implemented in

			**double Psi(const Eigen::Vector2d y)**
			**Eigen::Vector2d gradPsi(const Eigen::Vector2d y)**
			**double laplPsi(const Eigen::Vector2d y)**

		    Moreover, replace current signature of Jstar function as follows:

			**template <typename FUNCTOR>
			double Jstar(std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space,
            FUNCTOR &&u, const Eigen::Vector2d x)**
		   
			Adapt the description of input arguments. Note that I use a FUNCTOR for u, instead of the finite element vector u_N as input argument. I believe that I do not need to approximate u in the finite element space. Am I mistaken?

			Modify Hint accordingly!! There are wrong function signatures. Instead of Gval and Ggradval, add
			double G(Eigen::Vector2d x, Eigen::Vector2d y)

			and

			Eigen::Vector2d gradG(Eigen::Vector2d x, Eigen::Vector2d y)


- (3-11.g): replace the current function signature by

			**template <typename FUNCTOR>
			double stab_pointEval(
    			std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space,
    			FUNCTOR &&u, const Eigen::Vector2d x);**

			Note that I do not make use of a so called function Solve_LaplDir_LFE! I use directly the function u(x) given in the problem description. See also (3-11.f).


- (3-11.h): They shall adapt the main function, **stableevaluationatapoint_main.cc**.


PLEASE NOTE: The implementation for stableEval (Using the function JStar) does not converge. I could not figure out the error. I checked the formula for the laplacian of Psi several times, I believe it is correct. Open for discussion :)
