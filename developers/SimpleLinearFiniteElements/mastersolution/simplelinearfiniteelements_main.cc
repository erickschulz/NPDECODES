//#include <cstdlib>
#include <iostream>

#include "simplelinearfiniteelements.h"

/* SAM_LISTING_BEGIN_5 */
#define MESH "Square4"

int main()
{
  std::string meshfile = CURRENT_SOURCE_DIR "/../meshes/" MESH ".txt";

  SimpleLinearFiniteElements::TriaMesh2D square_mesh(meshfile);
  std::cout << "Mesh loaded " << std::endl;
  std::cout << "Mesh info: " << square_mesh.Coordinates.rows() << " vertices, "
            << square_mesh.Elements.rows() << " elements" << std::endl;

  // print both H1 and L2 errors and plot Mesh
  std::tuple<Eigen::VectorXd, double, double> solution = solve(square_mesh);

  std::cout << "L2-error:  " << std::get<1>(solution) << std::endl;
  std::cout << "H1s-error: " << std::get<2>(solution) << std::endl;

#if SOLUTION
  std::string meshplot = CURRENT_BINARY_DIR "/" MESH ".png";
  std::string meshfile_3d = CURRENT_BINARY_DIR "/" MESH "_3d.txt";
  std::string meshplot_3d = CURRENT_BINARY_DIR "/" MESH "_3d.png";

  // plot mesh
  std::system(("python3 -B " CURRENT_SOURCE_DIR "/plot_mesh.py " + meshfile + " " + meshplot).c_str());
  std::cout << "Generated " + meshplot << std::endl;

  // generate 3d mesh file from solution
  SimpleLinearFiniteElements::TriaMesh2D::addZComponent(meshfile, meshfile_3d, std::get<0>(solution));
  std::cout << "Generated " + meshfile_3d << std::endl;

  // plot the 3d mesh file
  std::system(("python3 -B " CURRENT_SOURCE_DIR "/plot_surf.py " + meshfile_3d + " " + meshplot_3d).c_str());
  std::cout << "Generated " + meshplot_3d << std::endl;
#else
  // To plot the mesh and your solution, uncomment the following:

  /*std::string meshplot = CURRENT_BINARY_DIR "/" MESH ".png";
  std::string meshfile_3d = CURRENT_BINARY_DIR "/" MESH "_3d.txt";
  std::string meshplot_3d = CURRENT_BINARY_DIR "/" MESH "_3d.png";

  // plot mesh
  std::system(("python3 -B " CURRENT_SOURCE_DIR "/plot_mesh.py " + meshfile + " " + meshplot).c_str());
  std::cout << "Generated " + meshplot << std::endl;

  // generate 3d mesh file from solution
  SimpleLinearFiniteElements::TriaMesh2D::addZComponent(meshfile, meshfile_3d, std::get<0>(solution));
  std::cout << "Generated " + meshfile_3d << std::endl;

  // plot the 3d mesh file
  std::system(("python3 -B " CURRENT_SOURCE_DIR "/plot_surf.py " + meshfile_3d + " " + meshplot_3d).c_str());
  std::cout << "Generated " + meshplot_3d << std::endl;*/
#endif

  return 0;
}
/* SAM_LISTING_END_5 */
