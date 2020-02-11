# NPDECODES
This repository contains the codes for the homework problems of the recurring course **Numerical Methods for Partial Differential Equations** at [ETH Zurich](https://ethz.ch/en.html). The course treats finite element methods (FEMs) using the C++ library [LehrFEM++](https://github.com/craffael/lehrfempp) and relies on the following material:
* [lecture notes](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE19.pdf)
* [homework problems](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/HOMEWORK/NPDEFL_Problems.pdf) 
Moreover, enrolled students can access the [moodle](https://moodle-app2.let.ethz.ch/course/view.php?id=12060) page of the course.

## Getting started
Open a terminal and type
```
git clone git@github.com:OliverRietmann/NPDECODES.git
cd NPDECODES/
mkdir build
cd build/
cmake ..
```
This will install LehrFEM++ and its dependencies into a folder `~/.hunter/`. To build a specific problem, say `BurgersEquation`, proceed as follows:
```
cd homeworks/BurgersEquation/
make
```
This will build from the source files in `NPDECODES/homeworks/BurgersEquation/`, where the subfolder `mysolution/` contains templates to be changed by the students. Recompilation is done by invoking `make` again. The following executables are generated:
* `./BurgersEquation_mastersolution`: Runs the mastersolution.
* `./BurgersEquation_test_mastersolution`: Runs unit tests on all important functions of the mastersolution.
* `./BurgersEquation_mysolution`: Runs the students code, i.e. the one in `mysolution/`.
* `./BurgersEquation_test_mysolution`: Runs unit tests the students code, i.e. the one in `mysolution/`.

Some remarks on building the codes:
* If you just clone the repository in this way, you can only work locally on your computer, since you have no permission to push to this remote repository. It may thus be useful to create your own copy (fork) of this repository. GitHub offers a [tutorial](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) on how to create such a fork.
* LehrFEM is already installed on the linux student computers in the ETH main building. To access this installation, you have to set the correct installation directory: In the terminal, run `export HUNTER_ROOT=/opt/libs/NumPDE` before running `cmake ..`. Then proceed as above. This sets the environment variable `HUNTER_ROOT` in the current terminal instance. If you start a new terminal, then you need to set it again. 
* There is two folders called `homeworks`. The source files are located in `NPDECODES/homeworks/<ProblemName>` and the binaries get created in `NPDECODES/build/homeworks/<ProblemName>`.
