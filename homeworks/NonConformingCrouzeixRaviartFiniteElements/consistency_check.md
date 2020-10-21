
THESE INCONSISTENCIES WERE FIXED:

- 2-14.h)  Tell them to actually implement lf::assemble::UniformFEDofHandler dof_handler() in the main file

- 2-14.q): replace **base::RefEl** with **lf::base::RefEl::**
		   replace **size_type** with **lf::assemble::size_type**
		   replace **dim_t** with **lf::assemble::dim_t**
		   replace **sub_idx_t** with **lf::base::sub_idx_t**
		   replace **scalar_type** with **double**
		   there is an additional comment before line 6 in EvalReferenceShapeFunctions()

- 2-14.r): replace **size_type** with **lf::assemble::size_type**
		   replace **scalar_type** with **double**
		   there is an additional comment before line 6 and before line 10 in GradientsReferenceShapeFunctions()

- 2-14.s): replace **scalar_type** with **double**

- 2-14.t): replace **scalar_type** with **double**
		   Is this subexercise really necesseary ?

- 2-14.u): replace **dof_handler.NoDofs()** with **dof_handler.NumDofs()**
		   replace **lf::uscalfe::MeshFunctionGlobal** with **lf::mesh::utils::MeshFunctionGlobal**
           replace **scalar_type** with **double**

- 2-14.v): replace **scalar_type** with **double**
		   replace **glb_idx_t** with **lf::assemble::glb_idx_t**

- 2-14.w): replace **size_type** with **lf::assemble::size_type**
		   in NUMPDE problems text file: const **lf::base::RandomAccessRange<const gdof_idx_t> cell_dof_idx()** ; whereas in problem code: **nonstd::span<const lf::assemble::gdof_idx_t>cell_dof_idx ()**

- I added in the codes all the comments **// TODO: task 2-14.xxx)** for every subexercise
