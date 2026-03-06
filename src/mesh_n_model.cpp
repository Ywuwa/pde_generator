#include "../headers/mesh_n_model.hpp"
#include "../headers/inout.hpp"
//=================================== INITIAL FUCNTIONS ===========================================
using abc = ABC_Flow;
double abc::initV1(const double x, const double y, const double z, const double t, const double Re)
{
  // u = (A*sin(kz) + C*cos(ky)) * exp(-tk/Re)  
  const double coef = std::exp((-t*abc::k) / Re);
  return ( abc::A * std::sin(abc::k * z) + abc::C * std::cos(abc::k * y) ) * coef;
}
double abc::initV2(const double x, const double y, const double z, const double t, const double Re)
{
  // v = (B*sin(kx) + A*cos(kz)) * exp(-tk/Re) 
  const double coef = std::exp((-t*abc::k) / Re);
  return ( abc::B * std::sin(abc::k * x) + abc::A * std::cos(abc::k * z) ) * coef;
}
double abc::initV3(const double x, const double y, const double z, const double t, const double Re)
{
  // w = (C*sin(ky) + B*cos(kx)) * exp(-tk/Re) 
  const double coef = std::exp((-t*abc::k) / Re);
  return ( abc::C * std::sin(abc::k * y) + abc::B * std::cos(abc::k * x) ) * coef;
}
double abc::initPress(
  const double x, const double y, const double z, const double t, const double Re)
{
  const double u = abc::initV1(x, y, z);
  const double v = abc::initV2(x, y, z);
  const double w = abc::initV3(x, y, z);
  // p = p0 - V^2 / 2 - P
  return abc::p0 - 0.5 * (u*u + v*v + w*w) - abc::P;
}
//=================================================================================================

//=================================== INITIAL CONDITIONS ==========================================
void initialConditions(
  std::vector<double>& feature, const uint initFuncInd, 
  const model_data& params, const double startTime)
{
  const uint domainPartition = params.domainPartition;
  const uint offsetY = domainPartition + 1;
  const uint offsetZ = (domainPartition+1) * (domainPartition+1);
  const double xStep (params.xLen / domainPartition);
  const double yStep (params.yLen / domainPartition);
  const double zStep (params.zLen / domainPartition);
  
  ABC_Flow functionSet;
  functionContainer<ABC_Flow> fC(&functionSet);
  fC.setFunctions();
  for (size_t k = 0; k < domainPartition + 1; k++)      // Z-Axis
  {
    for (size_t j = 0; j < domainPartition + 1; j++)    // Y-Axis
    {
      for (size_t i = 0; i < domainPartition + 1; i++)  // X-Axis
      {
        uint index (k*offsetZ + j*offsetY + i);
        feature[index] = (
          functionSet.*fC.indexedFunc[initFuncInd])
          (i*xStep, j*yStep, k*zStep, startTime, params.Reyn) ;
      }
    }
  }
}
//=================================================================================================

//==================================== PRECISE SOLUTION ===========================================
void compute_precise(
  const model_data& params,
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p)
{
  uint tick = 0;
  const auto dimSize ( params.domainPartition );
  const uint offsetY = dimSize + 1;
  const uint offsetZ = (dimSize+1) * (dimSize+1);

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize );
  const double hY (params.yLen / dimSize );
  const double hZ (params.zLen / dimSize );

  const std::string outputFuncFile = params.PATH;
  ABC_Flow functionSet;
  functionContainer<ABC_Flow> fC(&functionSet);
  fC.setFunctions();
  while (tick < params.timePartition + 1)
  {
    for (size_t k = 0; k < dimSize + 1; k++)      // Z-Axis
    {
      for (size_t j = 0; j < dimSize + 1; j++)    // Y-Axis
      {
        for (size_t i = 0; i < dimSize + 1; i++)  // X-Axis
        {
          uint index (k*offsetZ + j*offsetY + i);
          u[index] = (functionSet.*fC.indexedFunc[0])(i*hX, j*hY, k*hZ, tick*tau, params.Reyn);
          v[index] = (functionSet.*fC.indexedFunc[1])(i*hX, j*hY, k*hZ, tick*tau, params.Reyn);
          w[index] = (functionSet.*fC.indexedFunc[2])(i*hX, j*hY, k*hZ, tick*tau, params.Reyn);
          p[index] = (functionSet.*fC.indexedFunc[3])(i*hX, j*hY, k*hZ, tick*tau, params.Reyn);
        }
      }
    }
    funcOutput(outputFuncFile, "/helical_v1", std::to_string(tick), ".txt", u, params, false);
    funcOutput(outputFuncFile, "/helical_v2", std::to_string(tick), ".txt", v, params, false);
    funcOutput(outputFuncFile, "/helical_v3", std::to_string(tick), ".txt", w, params, false);
    funcOutput(outputFuncFile, "/helical_p", std::to_string(tick), ".txt", p, params, false);
    tick += 1;
  }
  std::cout << "final tick: " << tick << '\n';
}
//=================================================================================================