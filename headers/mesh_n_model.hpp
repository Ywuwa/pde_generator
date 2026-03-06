#include "settings.hpp"
#include <functional>
//=================================== INITIAL FUCNTIONS ===========================================
/*! \brief Function set consits of 3 functions for the velocity components and the pressure
 */
struct ABC_Flow
{
  double A = 1.0;
  double B = 0.1;
  double C = 0.1;
  double k = 3.0; // the vorticity of flow
  double P = 1.0; // the potential of volumetric forces
  double p0 = std::fabs(A) + std::fabs(B) + std::fabs(C) + P;
  double initV1(
    const double x, const double y, const double z, const double t = 0.0, const double Re = 100.0);
  double initV2(
    const double x, const double y, const double z, const double t = 0.0, const double Re = 100.0);
  double initV3(
    const double x, const double y, const double z, const double t = 0.0, const double Re = 100.0);
  double initPress(
    const double x, const double y, const double z, const double t = 0.0, const double Re = 100.0);  
};
/*! \brief Function Container accumulate 4 functions of chosen function Set
 */
template <typename fSet>
struct functionContainer : model_data
/*{
  std::vector<double (fSet::*)(
    const double, const double, const double, const double, const double )> indexedFunc;
  functionContainer()
  { 
    indexedFunc.emplace_back(&fSet::initV1);
    indexedFunc.emplace_back(&fSet::initV2);
    indexedFunc.emplace_back(&fSet::initV3);
    indexedFunc.emplace_back(&fSet::initPress);
  }
};*/
{
  std::vector<double (fSet::*)(
    const double, const double, const double, const double, const double )> indexedFunc;
  functionContainer(fSet* name) : instancePtr(name){ }
  void setFunctions()
  {
    indexedFunc.emplace_back(&fSet::initV1);
    indexedFunc.emplace_back(&fSet::initV2);
    indexedFunc.emplace_back(&fSet::initV3);
    indexedFunc.emplace_back(&fSet::initPress);
  }
  public:
    fSet *instancePtr; // pointer to an fSet instance
};
//=================================================================================================

//=================================== INITIAL CONDITIONS ==========================================
/*! \brief Applying initial conditions for functions
 *  \param[in] feature     - feature that needs to be initialized
 *  \param[in] initFuncInd - the required function index (0, 1, 2 or 3)
 *  \param[in] params      - the link to parameters structure
 *  \param[in] startTime   - initial condition time (default = 0.0)
 */
void initialConditions(
  std::vector<double>& feature, const uint initFuncInd, 
  const model_data& params, const double startTime = 0.0);
//=================================================================================================

//===================================== AXIS MESH INIT ============================================
/*! \brief axis mesh init function
 *  \param[in] axisMesh - mesh along <AXIS>
 *  \param[in] params   - the link to parameters structure
 *  \param[in] axisLen  - the length of the current axis
 */
void meshInit(std::vector<double>& axisMesh, const model_data& params, const double axisLen);
//=================================================================================================

//==================================== PRECISE SOLUTION ===========================================
/*! \brief Flow computation function
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute_precise(
  const model_data& params,
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================