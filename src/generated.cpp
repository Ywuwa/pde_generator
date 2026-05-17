
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
               const size_t dimSize,
               const model_data& params)
{
    const double Q1 = 1/params.Reyn;
  //std::cout << "explicit eq gen start" << '\n';
  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + i*offset_X);
      u1[index] = u[index] - tau*(((((((u[index] * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X))) + (u[index] * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((u[index] * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y))) + (v[index] * ((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y))))) + ((u[index] * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z))) + (w[index] * ((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z))))) + ((p[index + offset_X] + (-1 * p[index - offset_X])) * 1/(2*h_X))) + (-1 * (Q1 * (((((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X)) + (((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y))) + (((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z)))))));
    v1[index] = v[index] - tau*(((((((u[index] * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X))) + (v[index] * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((v[index] * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y))) + (v[index] * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y))))) + ((v[index] * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z))) + (w[index] * ((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z))))) + ((p[index + offset_Y] + (-1 * p[index - offset_Y])) * 1/(2*h_Y))) + (-1 * (Q1 * (((((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X)) + (((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y))) + (((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z)))))));
    w1[index] = w[index] - tau*(((((((u[index] * ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X))) + (w[index] * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((v[index] * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y))) + (w[index] * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y))))) + ((w[index] * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z))) + (w[index] * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z))))) + ((p[index + offset_Z] + (-1 * p[index - offset_Z])) * 1/(2*h_Z))) + (-1 * (Q1 * (((((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X)) + (((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y))) + (((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z)))))));
      }
    }
  }
  //std::cout << "explicit eq gen end" << '\n';
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
               const size_t dimSize,
               const model_data& params)
{
    const double Q1 = 1/params.Reyn;
  //std::cout << "implicit eq gen start" << '\n';
  for (size_t k = 2; k < dimSize-1; k++)      // Z-Axis
  {
    for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
    {
      for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
      {
      uint index (k*offset_Z + j*offset_Y + i*offset_X);
      const double indexVal_L00 = 1/(4*h_X*h_X);
            triplets0.emplace_back(index,index+(-2)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_L00);
    const double indexVal_0L0 = 1/(4*h_Y*h_Y);
            triplets0.emplace_back(index,index+(0)*offset_X+(-2)*offset_Y+(0)*offset_Z,indexVal_0L0);
    const double indexVal_00L = 1/(4*h_Z*h_Z);
            triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(-2)*offset_Z,indexVal_00L);
    const double indexVal_000 = (1/(4*h_Z*h_Z) * -1) + (1/(4*h_Z*h_Z) * -1) + (1/(4*h_Y*h_Y) * -1) + (1/(4*h_Y*h_Y) * -1) + (1/(4*h_X*h_X) * -1) + (1/(4*h_X*h_X) * -1);
            triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_000);
    const double indexVal_00R = 1/(4*h_Z*h_Z);
            triplets0.emplace_back(index,index+(0)*offset_X+(0)*offset_Y+(2)*offset_Z,indexVal_00R);
    const double indexVal_0R0 = 1/(4*h_Y*h_Y);
            triplets0.emplace_back(index,index+(0)*offset_X+(2)*offset_Y+(0)*offset_Z,indexVal_0R0);
    const double indexVal_R00 = 1/(4*h_X*h_X);
            triplets0.emplace_back(index,index+(2)*offset_X+(0)*offset_Y+(0)*offset_Z,indexVal_R00);
    
    B0[index] =  -(((((((((((((((((((u[index] * (((u[index + 2*offset_X] + (-1 * u[index])) + ((-1 * u[index]) + u[index - 2*offset_X])) * 1/(4*h_X*h_X))) + (((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)) * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + (u[index] * (((u[index + 2*offset_X] + (-1 * u[index])) + ((-1 * u[index]) + u[index - 2*offset_X])) * 1/(4*h_X*h_X)))) + (((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)) * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + (v[index] * (((v[index + 2*offset_Y] + (-1 * v[index])) + ((-1 * v[index]) + v[index - 2*offset_Y])) * 1/(4*h_Y*h_Y)))) + (((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)) * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + (v[index] * (((v[index + 2*offset_Y] + (-1 * v[index])) + ((-1 * v[index]) + v[index - 2*offset_Y])) * 1/(4*h_Y*h_Y)))) + (((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)) * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + (w[index] * (((w[index + 2*offset_Z] + (-1 * w[index])) + ((-1 * w[index]) + w[index - 2*offset_Z])) * 1/(4*h_Z*h_Z)))) + (((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) + (w[index] * (((w[index + 2*offset_Z] + (-1 * w[index])) + ((-1 * w[index]) + w[index - 2*offset_Z])) * 1/(4*h_Z*h_Z)))) + (((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)))) + (2 * (((u[index] * (((v[index + offset_X + offset_Y] + (-1 * v[index + offset_X - offset_Y])) + ((-1 * v[index - offset_X + offset_Y]) + v[index - offset_X - offset_Y])) * 1/(4*h_X*h_Y))) + (((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)) * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((v[index] * (((u[index + offset_X + offset_Y] + (-1 * u[index + offset_X - offset_Y])) + ((-1 * u[index - offset_X + offset_Y]) + u[index - offset_X - offset_Y])) * 1/(4*h_X*h_Y))) + (((u[index + offset_Y] + (-1 * u[index - offset_Y])) * 1/(2*h_Y)) * ((v[index + offset_X] + (-1 * v[index - offset_X])) * 1/(2*h_X))))))) + (2 * (((u[index] * (((w[index + offset_X + offset_Z] + (-1 * w[index + offset_X - offset_Z])) + ((-1 * w[index - offset_X + offset_Z]) + w[index - offset_X - offset_Z])) * 1/(4*h_X*h_Z))) + (((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) * ((u[index + offset_X] + (-1 * u[index - offset_X])) * 1/(2*h_X)))) + ((w[index] * (((u[index + offset_X + offset_Z] + (-1 * u[index + offset_X - offset_Z])) + ((-1 * u[index - offset_X + offset_Z]) + u[index - offset_X - offset_Z])) * 1/(4*h_X*h_Z))) + (((u[index + offset_Z] + (-1 * u[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_X] + (-1 * w[index - offset_X])) * 1/(2*h_X))))))) + (2 * (((v[index] * (((w[index + offset_Y + offset_Z] + (-1 * w[index + offset_Y - offset_Z])) + ((-1 * w[index - offset_Y + offset_Z]) + w[index - offset_Y - offset_Z])) * 1/(4*h_Y*h_Z))) + (((w[index + offset_Z] + (-1 * w[index - offset_Z])) * 1/(2*h_Z)) * ((v[index + offset_Y] + (-1 * v[index - offset_Y])) * 1/(2*h_Y)))) + ((w[index] * (((v[index + offset_Y + offset_Z] + (-1 * v[index + offset_Y - offset_Z])) + ((-1 * v[index - offset_Y + offset_Z]) + v[index - offset_Y - offset_Z])) * 1/(4*h_Y*h_Z))) + (((v[index + offset_Z] + (-1 * v[index - offset_Z])) * 1/(2*h_Z)) * ((w[index + offset_Y] + (-1 * w[index - offset_Y])) * 1/(2*h_Y))))))) + (-1 * (Q1 * (((((((((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X)))) + (((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X))) + (-1 * (((((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X)))) + (((u[index + offset_X] + (-2 * u[index])) + u[index - offset_X]) * 1/(h_X*h_X))))) * 1/(2*h_X*h_X*h_X)) + (((((((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y))) + (-1 * (((((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((u[index + offset_Y] + (-2 * u[index])) + u[index - offset_Y]) * 1/(h_Y*h_Y))))) * 1/(2*h_X*h_Y*h_Y))) + (((((((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z))) + (-1 * (((((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((u[index + offset_Z] + (-2 * u[index])) + u[index - offset_Z]) * 1/(h_Z*h_Z))))) * 1/(2*h_X*h_Z*h_Z)))))) + (-1 * (Q1 * (((((((((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X)))) + (((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X))) + (-1 * (((((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X)))) + (((v[index + offset_X] + (-2 * v[index])) + v[index - offset_X]) * 1/(h_X*h_X))))) * 1/(2*h_Y*h_X*h_X)) + (((((((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y))) + (-1 * (((((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((v[index + offset_Y] + (-2 * v[index])) + v[index - offset_Y]) * 1/(h_Y*h_Y))))) * 1/(2*h_Y*h_Y*h_Y))) + (((((((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z))) + (-1 * (((((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((v[index + offset_Z] + (-2 * v[index])) + v[index - offset_Z]) * 1/(h_Z*h_Z))))) * 1/(2*h_Y*h_Z*h_Z)))))) + (-1 * (Q1 * (((((((((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X)))) + (((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X))) + (-1 * (((((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X)) + (-2 * (((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X)))) + (((w[index + offset_X] + (-2 * w[index])) + w[index - offset_X]) * 1/(h_X*h_X))))) * 1/(2*h_Z*h_X*h_X)) + (((((((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y))) + (-1 * (((((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y)) + (-2 * (((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y)))) + (((w[index + offset_Y] + (-2 * w[index])) + w[index - offset_Y]) * 1/(h_Y*h_Y))))) * 1/(2*h_Z*h_Y*h_Y))) + (((((((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z))) + (-1 * (((((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z)) + (-2 * (((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z)))) + (((w[index + offset_Z] + (-2 * w[index])) + w[index - offset_Z]) * 1/(h_Z*h_Z))))) * 1/(2*h_Z*h_Z*h_Z)))))));
      }
    }
  }
  //std::cout << "implicit eq gen end" << '\n';
}
