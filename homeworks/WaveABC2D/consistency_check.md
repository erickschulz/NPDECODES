- (6-7.k): in class description on LaTeX file:
		   add **template <typename FUNC_RHO, typename FUNC_MU0, typename FUNC_NU0>** on top of class definition: that is

		   **template <typename FUNC_RHO, typename FUNC_MU0, typename FUNC_NU0>
			 class WaveABC2DTimestepper {
			 public:
			 ...
		   **
		   instead of having **template <typename FUNC_*>** inside class
		   
		   also remove comment **// The only public member function** 
		   add instead the additional member function
		   ** double energies();**

		   change in private members **double M_** to ** unsigned int M_** 

		   also remove any **FUNC_RHO &&rho** with **FUNC_RHO rho** 
		   (in constructor of WaveABC2DTIMESTEPPER)
		   
		   further below in the LaTeX file: replace **compGalerkinMat(fe_space_p, FUNC_ALPHA &&alpha, FUNC_GAMMA &&gamma, FUNC_BETA &&beta)** with **compGalerkinMat(fe_space_p, FUNC_ALPHA alpha, FUNC_GAMMA gamma, FUNC_BETA beta)**
		  
		  
		   also, where it reads "(I) The constructor": 
		   replace **template <typename FUNC_RHO>
		             WaveABC2DTimestepper::WaveABC2DTimestepper(
					 const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>
					 &fe_space_p,FUNC_RHO &&rho, unsigned i n t M, double T);** 

			with 
		   **template <typename FUNC_RHO, typename FUNC_MU0, typename FUNC_NU0>
		     WaveABC2DTimestepper<FUNC_RHO, FUNC_MU0, FUNC_NU0>::WaveABC2DTimestepper(
			   const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
			   FUNC_RHO rho, unsigned int M, double T)**
            
		
		   Moreover, where it reads "(II) The member function":
		   replace **template <typename FUNC_MU0>
		   			 Eigen::VectorXd WaveABC2DTimestepper::solveWaveABC2D(FUNC_MU0 &&u0);** 

		   with
		   **template <typename FUNC_RHO, typename FUNC_MU0, typename FUNC_NU0>
		     Eigen::VectorXd WaveABC2DTimestepper<FUNC_RHO, FUNC_MU0, FUNC_NU0>::solveWaveABC2D(
    		 FUNC_MU0 mu0, FUNC_NU0 nu0)**


  		    
		   Also remove description of **u0** and add description for **mu0** and **nu0**.

- (6-7.l)        Add that they should implement 

**template <typename FUNC_RHO, typename FUNC_MU0, typename FUNC_NU0>
double
    WaveABC2DTimestepper<FUNC_RHO, FUNC_MU0, FUNC_NU0>::energies()**;
    
    in the file waveabc2d.h.

I again assume that the codes will be updated accordingly in the solution files of LaTeX, for I added minor changes in the code.


