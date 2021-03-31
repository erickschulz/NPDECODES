#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <mgl2/mgl.h>
#include "meshgrid.hpp"

using namespace Eigen;

class Func {
    public:
    MatrixXd operator()( MatrixXd& X, MatrixXd Y, double gamma = (3+std::sqrt(3))/6 )
    {
        MatrixXcd Z(X.rows(), X.cols());
        Z.real() = X;
        Z.imag() = Y;

        // Define stability function
        MatrixXcd num = Z.array().square()*(0.5 - 2*gamma + gamma*gamma);
        num += MatrixXcd::Ones(Z.rows(),Z.cols()) + Z*(1 - 2*gamma);
        MatrixXcd den = (MatrixXcd::Ones(Z.rows(),Z.cols())-gamma*Z).array().square();
        MatrixXcd tmp = num.cwiseQuotient(den);

        MatrixXd S = tmp.array().abs();

        return S;
    }
};

int main()
{
    double gamma;
    std::cout << "Please enter a value for gamma: ";
    std::cin >> gamma;

    // Initialization
    VectorXd X_lin = VectorXd::LinSpaced(301, -15, +15);
    VectorXd Y_lin = VectorXd::LinSpaced(301, -15, +15);
    MatrixXd X, Y;
    meshgrid(X_lin, Y_lin, X, Y);

    Func S;
    MatrixXd S_ = S(X, Y, gamma);

    // Normalize results for plot
    S_ = (S_.array()/double(S_.maxCoeff())).matrix();

    mglData X_gl(X.rows(), X.cols(), X.data());
    mglData Y_gl(Y.rows(), Y.cols(), Y.data());
    mglData S_gl(S_.rows(), S_.cols(), S_.data());

    mglGraph gr;
    gr.SetRanges(X_gl.Minimal(), X_gl.Maximal(), Y_gl.Minimal(), Y_gl.Maximal());
    gr.Colorbar("kRryw");
    gr.Cont(S_gl);
    std::string title = "Stability Domain for 2-Stage SDIRK, gamma = " + std::to_string(gamma);
    gr.Title(title.c_str());
    gr.Axis(); gr.Label('x',"Re",0); gr.Label('y',"Im",0);
    gr.Tile(X_gl, Y_gl, S_gl, "kRryw");
    gr.WriteEPS("stabdomSDIRK");
}
