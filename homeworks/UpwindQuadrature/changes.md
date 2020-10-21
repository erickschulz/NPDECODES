- Created new helper function "opposite_velocity_directions", which determines the direction of -v(a^j) relative to the triangle.

REASON: Makes the remaining computation in the ElementMatrixProvider class a shorter and more readable. The separation also allows separate testing of this aspect of the implementation.

- Dropped implementation of separate assembly function

REASON: In Betl, data used during the local computations were passed to the local assembler from the assembly funcion using data objects, containing the parameters as fields. In lehrfempp this is no longer the case and parameters are usually passed to the ElementMatrixProvider via the constructor, thus the specialized assembly function is no longer required.

- Adapted representation of masses (vector -> MeshDataSet). The computation of the masses is now implemented in another helper function "initialize_masses".

REASON: During local assembly the ElementMatrixProvider has (as far as I know) no access to the index of the entity, thus an index based access to the masses is no longer possible. 
