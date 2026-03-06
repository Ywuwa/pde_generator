
#include "../headers/generated.hpp"
void generated(std::vector<double>& u,
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
