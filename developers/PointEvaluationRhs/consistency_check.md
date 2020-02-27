Problem 3-3: Dirichlet BVP with point-evaluation right-hand-side functional

- in PDF, 'functiona' should be 'functional' in title.

- in 3-3.d), e), f), g)(PDF): add **In the file pointevaluationrhs.cc**.

- in 3-3.g), add `If SOLUTION` tag.

- in 3-3.j)(PDF), In the file ~~pointEvaluation.cc~~ pointevaluationrhs.cc.

  And the code for definition of `class DeltaLocalVectorAssembler` is inconsistency with that in the `pointevaluationrhs.h`, it should be

  ```c++
  class DeltaLocalVectorAssembler{
  private:
    Eigen::Vector2d x_0;
    bool already_found;
  
  public:
    explicit DeltaLocalVectorAssembler(Eigen::Vector2d x)
        : x_0(x), already_found(false) {}
    bool isActive(const lf::mesh::Entity &entity) const{
      return (!already_found);
    }
    Eigen::VectorXd Eval(const lf::mesh::Entity &entity);
  };
  ```

- in 3-3.k), l): in the file ~~norms.cc~~ pointevaluationrhs_norms.cc

- in 3-3.m): in the file ~~PointEvaluation.cc~~ pointevaluationrhs.cc

- in `test/pointevaluationrhs_test.cc`: remove the cout.

