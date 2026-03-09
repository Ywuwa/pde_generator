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

//======================================== RESIDUALS ==============================================
/*! \brief Velocity residual computation function, based on FDA1 and FDA3
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double velocity_residual(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);
//=================================================================================================