# NPDECODES
This repository contains the codes for the homework problems of the recurring course **Numerical Methods for Partial Differential Equations** at [ETH Zurich](https://ethz.ch/en.html). The course treats finite element methods (FEMs) using the C++ library [LehrFEM++](https://github.com/craffael/lehrfempp) and relies on the following material:
* [lecture notes](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
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
This section is suited only for your own computer. To build the codes on the student computers of ETH see below. Open a terminal and type
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
    :   ├── mastersolution (folder containing source files)
        ├── mysolution     (folder containing source files, to be modified by you)
        ├── templates      (folder containing source files)
        :
```

## On the student computers
LehrFEM++ is already installed on the linux student computers in the ETH main building. To set up your local repository there, open a terminal an type:
```
cd /tmp
git clone git@github.com:erickschulz/NPDECODES.git
mv NPDECODES ~
cd ~/NPDECODES
mkdir build
cd build
export HUNTER_ROOT=/opt/libs/NumPDE
cmake ..
```
The first four lines are due to limited resources on the student computers. Setting the environment variable `export HUNTER_ROOT=/opt/libs/NumPDE` tells CMake where to look for the preinstalled libraries. This environment variable is local to your terminal, i.e. has to be redefined if you start a new terminal. Apart from this, you can use the folder `~/NPDECODES` in the same way you would for the approach in the previous section.
