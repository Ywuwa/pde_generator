
#include "../headers/generated.hpp"
void generated_time_eq(const uint offset_X,
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
      
      }
    }
  }
}

void generated_impl_eq(std::vector<double>& u,
               std::vector<Eigen::Triplet<double>>& triplets0,
               Eigen::VectorXd& B0,
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
      const double indexVal_L00 = ((1/h_X * -1) * Q1);
                triplets0.emplace_back(index,index+(-1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
    const double indexVal_000 = ((1/h_X * Q1) * -1) + (1/h_X * Q1) + (1/h_X * Q1) + 1/tau;
                triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
    
    B0[index] = 1/tau*u[index] -(0);
      }
    }
  }
}
