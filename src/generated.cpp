
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
               std::vector<double>& v,
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
  

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      const double indexVal_L00 = (1/(2*h_X) * -1);
                triplets0.emplace_back(index,index+(-1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
    const double indexVal_000 = 1/tau;
                triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
    const double indexVal_R00 = 1/(2*h_X);
                triplets0.emplace_back(index,index+(1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_R00);
    
    B0[index] = 1/tau*u[index] -((((v[index + offset_Z + offset_Z] + (-1 * v[index + offset_Z - offset_Z])) + ((-1 * v[index - offset_Z + offset_Z]) + v[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)));
      }
    }
  }
}
