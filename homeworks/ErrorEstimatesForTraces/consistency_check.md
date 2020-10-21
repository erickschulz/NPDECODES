- 3.5-i) replace **BOOST_ASSERT()** with **assert()**
		 
		 replace **auto edges_predicate = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
    	 return bd_flags(edge);** with **auto edges_predicate = [&bd_flags](const lf::mesh::Entity *edge) -> bool {
    	 return bd_flags(*edge);**
		
		 replace **for (const lf::mesh::Entity &edge : mesh_p->Entities(1))** with **for (const lf::mesh::Entity *edge : mesh_p->Entities(1))**
		 
		 replace **if (bd_flags(edge))** with **if (bd_flags(*edge))**

- 3.5-j) replace **boost::filesystem::path here = __FILE__;
				   std::string filename = "/meshes/hex" + idx_str + ".msh";
				   auto mesh_path = here.parent_path().parent_path() / filename;**
	     with **std::string mesh_file = CURRENT_SOURCE_DIR "/../meshes/hex" + std::to_string(i) + ".msh";**

		 and replace **const lf::io::GmshReader reader(std::move(mesh_factory),mesh_path.string());** with const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);**

		 replace **dofh.NoDofs()** with **dofh.NumDofs()**

		 there is an **std::abs()** missing in code line 25 in lecture document: use
		 double error = std::abs(bd_functional_val - 2.081541059732923);

		 replace in code line 17 in lecture documents **results(i-1, 1) = N_dofs** with **results(i-1, 0) = N_dofs**
		 replace in code line 26 in lecture document **results(i-1, 0) = error** with **resuls(i-1, 1) = error**

- Naming convention: replace filenames in description as follows: all occurences of 
					**tee_lapl_robin_assembly.* ** shall be replaced with **teelaplrobinassembly.* **
					**trace_error_estimate_main.cc** shall be replaced with **errorestimatesfortraces_main.cc**

							

		 
