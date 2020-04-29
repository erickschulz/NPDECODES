- 3-8.c) rename filenames **OutputImpedanceBVP.* **to **outputimpedancebvp.* **
		 
		 replace code line 8 in problem solution:  **Eigen::VectorXd discrete_solution;** with Eigen::VectorXd discrete_solution(N_dofs);**
		 
		 replace **coord_t** with ** Eigen::Vector2d **
		 
	 	 replace  **lf::assemble::fix_flagged_solution_comp_alt-<double>()** with ** lf::assemble::FixFlaggedSolutionCompAlt<double>()**
		 
- 3-8.h) rename filenames **OutputImpedanceBVP.cc** to **outputimpedancebvp.cc**		
 		 replace code line 34 in problem solution: **for (const lf::mesh::Entity &edge : mesh_p->Entities(1))** with **for (const lf::mesh::Entity *edge : mesh_p->Entities(1))**

		 replace **if (edges_predicate_RobinBC(edge))** with **if (edges_predicate_RobinBC(*edge))**

		 replace ** auto dof_idx = dofh.GlobalDofIndices(edge);** with ** auto dof_idx = dofh.GlobalDofIndices(*edge);**
		
		 replace ** assert(dofh.NumLocalDofs(edge) == 2);** with **assert(dofh.NumLocalDofs(*edge) == 2);**


