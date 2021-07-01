#ifndef SAMPLE_MESHFUNCTION_H
#define SAMPLE_MESHFUNCTION_H

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <cmath>
#include <memory>
#include <fstream>

namespace ConvectionDiffusion{

template<typename MF>
double EvaluateMeshFunction(std::shared_ptr<const lf::mesh::Mesh> mesh_p, MF mf, Eigen::Vector2d global, double tol=std::pow(10,-10)){
    for(const lf::mesh::Entity* entity_p: mesh_p->Entities(0)){
        //compute geometric information about the cell
        const lf::geometry::Geometry* geo_p = entity_p->Geometry();
        Eigen::MatrixXd corners = lf::geometry::Corners(*geo_p);
        
        //transform global coordinates to local coordinates on the cell
        Eigen::Matrix2d A;
        A << corners.col(1) - corners.col(0), corners.col(2) - corners.col(0);
        Eigen::Vector2d b;
        b << global - corners.col(0);
        Eigen::Vector2d loc = A.fullPivLu().solve(b); 

        //evaluate meshfunction, if local coordinates lie in the reference triangle
        if(loc(0) >= 0 - tol && loc(1) >= 0 - tol && loc(0) + loc(1) <= 1 + tol){
           return mf(*entity_p,loc)[0];
        }

    }
    //TODO: Add assert here
    return 0.0;
}

template<typename CURVE, typename MF>
void SampleMeshFunction(std::string file_name,std::shared_ptr<const lf::mesh::Mesh> mesh_p, CURVE gamma, MF mf, int N){
    //open file
    std::ofstream file;
    file.open(file_name); 

    //time points to sample along gamma
    Eigen::VectorXd sample_times = Eigen::VectorXd::LinSpaced(N,0.0,1.0);


    for(int i = 0; i < N; ++i){
        double t = sample_times(i);
        double eval = EvaluateMeshFunction(mesh_p,mf,gamma(t));
        file << t << ", " << eval << "\n";
    }  
    file.close();
}


} //namespace ConvectionDiffusion



#endif //SAMPLE_MESHFUNCTION_H