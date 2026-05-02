
#include "../headers/generated.hpp"
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
               const size_t dimSize)
{
    const double Q1 = 1/Re;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      u1[index] = u[index] + tau*((((((u[index] * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X))) + (v[index] * ((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)))) + (w[index] * ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)))) + ((p[index + offset_X] + (-1 * p[index - offset_X])) * 1/(2*h_X))) + (-1 * (Q1 * (((((u[index + offset_X + offset_X] + (-1 * u[index + offset_X - offset_X])) + ((-1 * u[index - offset_X + offset_X]) + u[index - offset_X - offset_X])) * 1/(4*h_X*h_X)) + (((u[index + offset_Y + offset_Y] + (-1 * u[index + offset_Y - offset_Y])) + ((-1 * u[index - offset_Y + offset_Y]) + u[index - offset_Y - offset_Y])) * 1/(4*h_Y*h_Y))) + (((u[index + offset_Z + offset_Z] + (-1 * u[index + offset_Z - offset_Z])) + ((-1 * u[index - offset_Z + offset_Z]) + u[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)))))));
    v1[index] = v[index] + tau*((((((u[index] * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X))) + (v[index] * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + (w[index] * ((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)))) + ((p[index + offset_Y] + (-1 * p[index - offset_Y])) * 1/(2*h_Y))) + (-1 * (Q1 * (((((v[index + offset_X + offset_X] + (-1 * v[index + offset_X - offset_X])) + ((-1 * v[index - offset_X + offset_X]) + v[index - offset_X - offset_X])) * 1/(4*h_X*h_X)) + (((v[index + offset_Y + offset_Y] + (-1 * v[index + offset_Y - offset_Y])) + ((-1 * v[index - offset_Y + offset_Y]) + v[index - offset_Y - offset_Y])) * 1/(4*h_Y*h_Y))) + (((v[index + offset_Z + offset_Z] + (-1 * v[index + offset_Z - offset_Z])) + ((-1 * v[index - offset_Z + offset_Z]) + v[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)))))));
    w1[index] = w[index] + tau*((((((u[index] * ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X))) + (v[index] * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y)))) + (w[index] * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) + ((p[index + offset_Z] + (-1 * p[index - offset_Z])) * 1/(2*h_Z))) + (-1 * (Q1 * (((((w[index + offset_X + offset_X] + (-1 * w[index + offset_X - offset_X])) + ((-1 * w[index - offset_X + offset_X]) + w[index - offset_X - offset_X])) * 1/(4*h_X*h_X)) + (((w[index + offset_Y + offset_Y] + (-1 * w[index + offset_Y - offset_Y])) + ((-1 * w[index - offset_Y + offset_Y]) + w[index - offset_Y - offset_Y])) * 1/(4*h_Y*h_Y))) + (((w[index + offset_Z + offset_Z] + (-1 * w[index + offset_Z - offset_Z])) + ((-1 * w[index - offset_Z + offset_Z]) + w[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)))))));
      }
    }
  }
}

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
               const size_t dimSize)
{
    const double Q1 = 1/Re;

  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + offset_X);
      const double indexVal_L00 = 1/(4*h_X*h_X) + (1/(4*h_X*h_X) * -1) + (1/(4*h_X*h_X) * -1);
            triplets0.emplace_back(index,index+(-1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
    const double indexVal_0L0 = 1/(4*h_Y*h_Y) + (1/(4*h_Y*h_Y) * -1) + (1/(4*h_Y*h_Y) * -1);
            triplets0.emplace_back(index,index+(0)*offset_X+(-1)*offset_Y+(0)*offset_Z,indexVal_0L0);
    const double indexVal_00L = 1/(4*h_Z*h_Z) + (1/(4*h_Z*h_Z) * -1) + (1/(4*h_Z*h_Z) * -1);
            triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(-1)*offset_Z,indexVal_00L);
    const double indexVal_00R = 1/(4*h_Z*h_Z);
            triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(1)*offset_Z,indexVal_00R);
    const double indexVal_0R0 = 1/(4*h_Y*h_Y);
            triplets0.emplace_back(index,index+(0)*offset_X+(1)*offset_Y+(0)*offset_Z,indexVal_0R0);
    const double indexVal_R00 = 1/(4*h_X*h_X);
            triplets0.emplace_back(index,index+(1)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_R00);
    
    B0[index] =  -((((((((((((((((u[index] * (((u[index + offset_X + offset_X] + (-1 * u[index + offset_X - offset_X])) + ((-1 * u[index - offset_X + offset_X]) + u[index - offset_X - offset_X])) * 1/(4*h_X*h_X))) + (((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)) * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + (v[index] * (((v[index + offset_Y + offset_Y] + (-1 * v[index + offset_Y - offset_Y])) + ((-1 * v[index - offset_Y + offset_Y]) + v[index - offset_Y - offset_Y])) * 1/(4*h_Y*h_Y)))) + (((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)) * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + (w[index] * (((w[index + offset_Z + offset_Z] + (-1 * w[index + offset_Z - offset_Z])) + ((-1 * w[index - offset_Z + offset_Z]) + w[index - offset_Z - offset_Z])) * 1/(4*h_Z*h_Z)))) + (((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) + (v[index] * (((u[index + offset_X + offset_Y] + (-1 * u[index + offset_X - offset_Y])) + ((-1 * u[index - offset_X + offset_Y]) + u[index - offset_X - offset_Y])) * 1/(4*h_X*h_Y)))) + (((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)) * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X)))) + (u[index] * (((w[index + offset_Z + offset_X] + (-1 * w[index + offset_Z - offset_X])) + ((-1 * w[index - offset_Z + offset_X]) + w[index - offset_Z - offset_X])) * 1/(4*h_Z*h_X)))) + (((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X)) * ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)))) + (w[index] * (((v[index + offset_Y + offset_Z] + (-1 * v[index + offset_Y - offset_Z])) + ((-1 * v[index - offset_Y + offset_Z]) + v[index - offset_Y - offset_Z])) * 1/(4*h_Y*h_Z)))) + (((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y)))) + (-1 * (Q1 * (((((((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)) + (-1 * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((-1 * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X))) + ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) * 1/(4*h_X*h_X)) + (((((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)) + (-1 * ((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)))) + ((-1 * ((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y))) + ((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)))) * 1/(4*h_X*h_Y))) + (((((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)) + (-1 * ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)))) + ((-1 * ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z))) + ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)))) * 1/(4*h_X*h_Z)))))) + (-1 * (Q1 * (((((((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X)) + (-1 * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X)))) + ((-1 * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X))) + ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X)))) * 1/(4*h_Y*h_X)) + (((((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)) + (-1 * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + ((-1 * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y))) + ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) * 1/(4*h_Y*h_Y))) + (((((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)) + (-1 * ((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)))) + ((-1 * ((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z))) + ((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)))) * 1/(4*h_Y*h_Z)))))) + (-1 * (Q1 * (((((((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X)) + (-1 * ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X)))) + ((-1 * ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X))) + ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X)))) * 1/(4*h_Z*h_X)) + (((((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y)) + (-1 * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y)))) + ((-1 * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y))) + ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y)))) * 1/(4*h_Z*h_Y))) + (((((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) + (-1 * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) + ((-1 * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z))) + ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) * 1/(4*h_Z*h_Z)))))));
      }
    }
  }
}
