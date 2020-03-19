- 6-5.g) in the problem description: replace 


<----------
**Eigen::SparseMatrix<double> assembleGalerkinMatrix(
	std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p, 
	FUNC_ALPHA&& alpha, FUNC_GAMMA&& gamma, FUNC_BETA&& beta)** with 

**Eigen::SparseMatrix<double> assembleGalerkinMatrix(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_p,
    FUNC_ALPHA alpha, FUNC_GAMMA gamma, FUNC_BETA beta)**
------> decided to stick with && and fixed code instead
		note: corresponing code needs adaptions (that is, replace all **lf::uscalfe::MeshFunctionGlobal()** with ** lf::mesh::utils::MeshFunctionGlobal()**

- 6-5.h) in problem description: place **template <typename FUNCTION>** outside class definition

- 6-5.l) code needs adaptions (especially mesh_factory (not via boost))

