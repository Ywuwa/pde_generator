
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
      const double C1 = 5;
        const double C2 = 8;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      u1[index] = u[index] + tau*((C1 * (u[index + offset_X] - u[index - offset_X]) / (2*h_X)));
    v1[index] = v[index] + tau*((C2 * (v[index + offset_Y] - v[index - offset_Y]) / (2*h_Y)));
      }
    }
  }
}

void generated_impl_eq(std::vector<double>& p,
               std::vector<double>& u,
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
      const double C1 = 5;
        const double C2 = 8;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      const double indexVal_00L = ((1/(2*h_Z) * -1) * u[index]);
                triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(-1)*offset_Z,indexVal_00L);
    const double indexVal_000 = ((1/(2*h_Z) * u[index - offset_Z]) * -1) + (1/(2*h_Z) * u[index + offset_Z]);
                triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
    const double indexVal_00R = (1/(2*h_Z) * u[index]);
                triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(1)*offset_Z,indexVal_00R);
    
    B0[index] = -(((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)));
      }
    }
  }
}
