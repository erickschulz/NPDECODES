# A boundary value problem modelling stationary heat conduction / Unstable BVP

Here we solve the heat equation for a triangular domain with a non-continuous boundary condition,
depending on where the domain is. We check convergence properties of the solutions for the possible
cases and find differences.

The file `h1.txt` contains the output of the executable, which is the values of the H1 seminorm of
the solutions for different refinement levels for different locations of the domain, as well as the 
difference of the seminorms to the most refined level.

The script `plot.py` visualises the output and the results are shown in `h1_seminorms(differences).png`.

