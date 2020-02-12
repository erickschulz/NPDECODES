# NPDECODES
This repository contains the codes for the homework problems of the recurring course **Numerical Methods for Partial Differential Equations** at [ETH Zurich](https://ethz.ch/en.html). The course treats finite element methods (FEMs) using the C++ library [LehrFEM++](https://github.com/craffael/lehrfempp) and relies on the following material:
* [lecture notes](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE19.pdf)
* [homework problems](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/HOMEWORK/NPDEFL_Problems.pdf)

Moreover, enrolled students can access the [moodle](https://moodle-app2.let.ethz.ch/course/view.php?id=12060) page of the course.

## Requirements
Currently, only UNIX based operating systems are supported. Moreover, you need to have the following installed on your machine:
* C++17 compiler (e.g. gcc, clang)
* CMake (at least VERSION 3.10)
* python3
* A reader for .vtk files (e.g. paraview)
* git (not strictly needed, you could also download the repo as .zip file)

## Getting started
Open a terminal and type
```
git clone git@github.com:erickschulz/NPDECODES.git
cd NPDECODES/
mkdir build
cd build/
cmake ..
```
This will install LehrFEM++ and its dependencies into a folder `~/.hunter/`. To build a specific problem, say `TestQuadratureRules`, proceed as follows:
```
cd homeworks/TestQuadratureRules/
make
```
This will build from the source files in `NPDECODES/homeworks/TestQuadratureRules/`, where the subfolder `mysolution/` contains templates to be changed by the students. Recompilation is done by invoking `make` again. The following executables are generated:
* `./TestQuadratureRules_mastersolution`: Runs the mastersolution.
* `./TestQuadratureRules_test_mastersolution`: Runs unit tests on all important functions of the mastersolution.
* `./TestQuadratureRules_mysolution`: Runs the students code, i.e. the one in `mysolution/`.
* `./TestQuadratureRules_test_mysolution`: Runs unit tests the students code, i.e. the one in `mysolution/`.

There is two folders called `homeworks/`. One contains the source files and one contains the executables:
```
.
├── build (was created by you)
│   ├── homeworks
│   :   ├── TestQuadratureRules
│       :   ├── TestQuadratureRules_mastersolution      (executable)
│           ├── TestQuadratureRules_mysolution          (executable)
│           ├── TestQuadratureRules_test_mastersolution (executable)
│           ├── TestQuadratureRules_test_mysolution     (executable)
│           :
│
├── homeworks
:   ├── TestQuadratureRules
    :   ├── mastersolution (contains source files)
        ├── mysolution     (contains source files, to be modified by you)
        ├── templates      (contains source files)
        :
```

## Further remarks
* If you just clone the repository in this way, you can only work locally on your computer, since you have no permission to push to this remote repository. It may thus be useful to create your own copy (fork) of this repository. GitHub offers a [tutorial](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) on how to create such a fork.
* LehrFEM++ is already installed on the linux student computers in the ETH main building. To access this installation, you have to set the correct installation directory: In the terminal, run `export HUNTER_ROOT=/opt/libs/NumPDE` before running `cmake ..`. Then proceed as above. This sets the environment variable `HUNTER_ROOT` in the current terminal instance. If you start a new terminal, then you need to set it again. 
