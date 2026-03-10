#include "settings.hpp"

//============================== FLOW COMPUTATION via FDA1/FDA3 ===================================
/*! \brief Flow computation function in the cubic domain, based on FDA1 or FDA3
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute_cube(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================