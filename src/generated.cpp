
#include "../headers/generated.hpp"
void generated_time_eq(std::vector<double>& u,
               std::vector<double>& u1,
               std::vector<double>& v,
               std::vector<double>& v1,
               const uint offset_X,
               const uint offset_Y,
               const uint offset_Z,
               const double h_X,
               const double h_Y,
               const double h_Z,
               const double tau,
               const size_t dimSize)
{
const double Q1 = 5;
const double Q2 = 6;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      u1[index] = u[index] + tau*( - (((1/(2*h_X) * v) * u[index + 1*offset_X]) + ((-1/(2*h_X) * v) * u[index + -1*offset_X]) + ((u * 1/(2*h_X)) * v[index + 1*offset_X]) + ((u * -1/(2*h_X)) * v[index + -1*offset_X])));
      }
    }
  }
}

void generated_impl_eq(std::vector<double>& u,
               std::vector<Eigen::Triplet<double>>& triplets0,
               Eigen::VectorXd& B0,
               std::vector<double>& u,
               std::vector<double>& v,
               std::vector<Eigen::Triplet<double>>& triplets1,
               Eigen::VectorXd& B1,
               const uint offset_X,
               const uint offset_Y,
               const uint offset_Z,
               const double h_X,
               const double h_Y,
               const double h_Z,
               const double tau,
               const size_t dimSize)
{
const double Q1 = 5;
const double Q2 = 6;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      const double indexVal_000 = (1/tau + (Q1 * (1/h_X * -1/h_Y)));
triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
const double indexVal_0R0 = (Q1 * (1/h_X * 1/h_Y));
triplets0.emplace_back(index,index+(0)*offset_X+(1)*offset_Y+(0)*offset_Z,indexVal_0R0);
const double indexVal_LR0 = (Q1 * (-1/h_X * 1/h_Y));
triplets0.emplace_back(index,index+(-1)*offset_X+(1)*offset_Y+(0)*offset_Z,indexVal_LR0);
const double indexVal_L00 = (Q1 * (-1/h_X * -1/h_Y));
triplets0.emplace_back(index,index+(-1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
B0[index] =  - ((Q1 + (2 * Q2)));
const double indexVal_00R = 1/h_Z;
triplets1.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(1)*offset_Z,indexVal_00R);
const double indexVal_000 = -1/h_Z;
triplets1.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
B1[index] =  - ((1/h_Z * v[index + 1*offset_Z]) + (-1/h_Z * v[index]));
      }
    }
  }
}
