
#include "settings.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
void generated_time_eq(std::vector<double>& p,
               std::vector<double>& u,
               std::vector<double>& v,
               std::vector<double>& w,
               std::vector<double>& u1,
               std::vector<double>& v1,
               std::vector<double>& w1,
               const uint offset_X,
               const uint offset_Y,
               const uint offset_Z,
               const double h_X,
               const double h_Y,
               const double h_Z,
               const double tau,
               const size_t dimSize,
               const model_data& params);
void generated_impl_eq(std::vector<double>& u,
               std::vector<double>& v,
               std::vector<double>& w,
               std::vector<Eigen::Triplet<double>>& triplets0,
               Eigen::VectorXd& B0,
               const uint offset_X,
               const uint offset_Y,
               const uint offset_Z,
               const double h_X,
               const double h_Y,
               const double h_Z,
               const double tau,
               const size_t dimSize,
               const model_data& params);
