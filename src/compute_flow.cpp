 
#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/compute_flow.hpp"
#include "../headers/generated.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

//======================================== RESIDUALS ==============================================

/*
<var>   - exact value at knot
<var>1  - estimated solution at knot
*/
double velocity_residual(const uint offset_X,
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
  double vectorResidual (0.0);
  //---------------------------- inner knots --------------------------------
  for (size_t k = 1; k < dimSize; k++)      // Z-Axis
  {
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index (k*offset_Z + j*offset_Y + i);
        double resTerm (0.0); // residual term
        //! Insert precise values to the scheme, take the difference with the estimated values
        
      }
    }
  }
  //-------------------------------------------------------------------------
  return std::sqrt(vectorResidual);
}

//=================================================================================================

//======================================= FLOW COMPUTATION ========================================
void compute_cube(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p0)
{
  uint tick = 1; // number of time step
  const auto dimSize ( params.domainPartition );  // 1-dimension size
  const uint offsetX = 1;
  const uint offsetY = dimSize + 1;               // j+1 component (i+1 component is just [index+1])
  const uint offsetZ = (dimSize+1) * (dimSize+1); // k+1 component

  const double tau ( params.duration / params.timePartition ); // time step
  const double hX (params.xLen / dimSize ); // x-step
  const double hY (params.yLen / dimSize ); // y-step
  const double hZ (params.zLen / dimSize ); // z-step

  std::ofstream outputFile(params.PATH + params.PATH_log, std::ios::out);
  std::ofstream outputResidualFile(params.PATH + params.PATH_residual, std::ios::out);
  const std::string outputFuncFile = params.PATH;
  // vector size
  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize, 0.0);
  std::vector<double> v1(vecSize, 0.0);
  std::vector<double> w1(vecSize, 0.0);
  Eigen::VectorXd p(vecSize);
  for (size_t i = 0; i < vecSize; i++) p[i] = p0[i];

  while (tick < params.timePartition + 1)
  {
    //! velocity exact
    std::vector<double> uExac(vecSize, 0.0);
    initialConditions(uExac, 0, params, tick*tau);
    std::vector<double> vExac(vecSize, 0.0);
    initialConditions(vExac, 1, params, tick*tau);
    std::vector<double> wExac(vecSize, 0.0);
    initialConditions(wExac, 2, params, tick*tau);
    std::vector<double> pExac(vecSize, 0.0);
    initialConditions(pExac, 3, params, tick*tau);

    //! velocity compute
    //---------------------------- inner knots --------------------------------
    generated_time_eq(offsetX, offsetY, offsetZ, hX, hY, hZ, tau, dimSize);
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0 / Z = MAX; Neiman's condition dw/dz = 0
    for (size_t j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];

        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
      }
    }
    //! XZ-plane Y = 0 / Y = MAX; Neiman's condition dv/dy = 0
    for (size_t k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (size_t i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];

        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];
      }
    }
    //! YZ-plane X = 0 / X = MAX; Neiman's condition du/dx = 0
    for (size_t k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (size_t j = 1; j < dimSize; j++)  // X-Axis
      {
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        w1[index1] = wExac[index1];
        w1[index2] = wExac[index2];
        v1[index1] = vExac[index1];
        v1[index2] = vExac[index2];

        u1[index1] = uExac[index1];
        u1[index2] = uExac[index2];
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0 || MAX); dv/dy = 0; dw/dz = 0
    for (size_t i = 1; i < dimSize; i++)
    {
      // Y = Z = 0
      uint index1 (i);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];
    }
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (size_t j = 1; j < dimSize; j++)
    {
      // X = Z = 0
      uint index1 (j*offsetY);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];
    }
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (size_t k = 1; k < dimSize; k++)
    {
      // X = Y = 0
      uint index1 (k*offsetZ);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);

      w1[index1] = wExac[index1];
      w1[index2] = wExac[index2];
      w1[index3] = wExac[index3];
      w1[index4] = wExac[index4];

      u1[index1] = uExac[index1];
      u1[index2] = uExac[index2];
      u1[index3] = uExac[index3];
      u1[index4] = uExac[index4];

      v1[index1] = vExac[index1];
      v1[index2] = vExac[index2];
      v1[index3] = vExac[index3];
      v1[index4] = vExac[index4];
    }
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 0
    uint index (0);
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    u1[index] = uExac[index];
    v1[index] = vExac[index];
    w1[index] = wExac[index];

    //-------------------------------------------------------------------------
    // move data from upper time layer
    u.clear(); u = std::move(u1);
    v.clear(); v = std::move(v1);
    w.clear(); w = std::move(w1);
    u1.resize(vecSize); v1.resize(vecSize); w1.resize(vecSize);
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //! initialization of equation entities (Ax = b)
    Eigen::SparseMatrix<double> A(vecSize, vecSize);
    Eigen::VectorXd B0(vecSize);
    std::vector<Eigen::Triplet<double>> triplets0; // entities for filling a sparse matrix

    //! pressure compute
    //---------------------------- inner knots --------------------------------
    generated_impl_eq(u0, triplets0, B0, offsetX, offsetY, offsetZ, hX, hY, hZ, tau, dimSize);
    //-------------------------------------------------------------------------

    //---------------------------- border knots -------------------------------
    //! XY-plane Z = 0,1 / Z = MAX-1,MAX; Neiman's condition dp/dz = 0
    for (size_t j = 2; j < dimSize-1; j++)    // Y-Axis
    {
      for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
      {
        uint index1 (j*offsetY + i);
        uint index2 (dimSize*offsetZ + j*offsetY + i);
        uint index3 (offsetZ + j*offsetY + i);
        uint index4 ((dimSize-1)*offsetZ + j*offsetY + i);
        // matrix construct via EIGEN triplets 
        // behind low border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }
    }
    //! XZ-plane Y = 0,1 / Y = MAX-1,MAX; Neiman's condition dp/dy = 0
    for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
    {
      for (size_t i = 2; i < dimSize-1; i++)  // X-Axis
      {
        uint index1 (k*offsetZ + i);
        uint index2 (k*offsetZ + dimSize*offsetY + i);
        uint index3 (k*offsetZ + offsetY + 1);
        uint index4 (k*offsetZ + (dimSize-1)*offsetY + i);
        // matrix construct via EIGEN triplets 
        // behind Y0 border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }
    }
    
    //! YZ-plane X = 0,1 / X = MAX-1,MAX; Neiman's condition dp/dx = 0
    for (size_t k = 2; k < dimSize-1; k++)    // Z-Axis
    {
      for (size_t j = 2; j < dimSize-1; j++)  // X-Axis
      {
        uint index1 (k*offsetZ + j*offsetY);
        uint index2 (k*offsetZ + j*offsetY + dimSize);
        uint index3 (k*offsetZ + j*offsetY + 1);
        uint index4 (k*offsetZ + j*offsetY + dimSize-1);
        // matrix construct via EIGEN triplets (flow out)
        // behind X0 border
        triplets0.emplace_back(index3,index3, 1.0);
        // low border
        triplets0.emplace_back(index1,index1, 1.0);
        // behind up border
        triplets0.emplace_back(index4,index4, 1.0);
        // up border
        triplets0.emplace_back(index2,index2, 1.0);
        B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
      }
    }
    //-------------------------------------------------------------------------

    //---------------------------- edge knots ---------------------------------
    //! Y = Z = (0,1 || MAX,MAX-1); 
    for (size_t i = 2; i < dimSize-1; i++)
    {
      // Y = Z = 0
      uint index1 (i);
      triplets0.emplace_back(index1,index1, 1.0);
      // Y = MAX, Z = 0
      uint index2 (dimSize*offsetY + i);
      triplets0.emplace_back(index2,index2, 1.0);
      // Y = 0, Z = MAX
      uint index3 (dimSize*offsetZ + i);
      triplets0.emplace_back(index3,index3, 1.0);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + dimSize*offsetY + i);
      triplets0.emplace_back(index4,index4, 1.0);

      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // Y = Z = 1
      index1 = offsetZ + offsetY + i;
      triplets0.emplace_back(index1,index1, 1.0);
      // Y = MAX-1, Z = 1
      index2 = (dimSize-1)*offsetY + offsetZ + i;
      triplets0.emplace_back(index2,index2, 1.0);
      // Y = 1, Z = MAX-1
      index3 = (dimSize-1)*offsetZ + offsetY + i;
      triplets0.emplace_back(index3,index3, 1.0);
      // Y = MAX-1, Z = MAX-1
      index4 = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + i;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }
    //! X = Z = (0 || MAX); du/dx = 0, dw/dz = 0
    for (size_t j = 2; j < dimSize-1; j++)
    {
      // X = Z = 0
      uint index1 (j*offsetY);
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX, Z = 0
      uint index2 (j*offsetY + dimSize);
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 0, Z = MAX
      uint index3 (dimSize*offsetZ + j*offsetY);
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX, Z = MAX
      uint index4 (dimSize*offsetZ + j*offsetY + dimSize);
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // X = Z = 1
      index1 = offsetZ + j*offsetY + 1;
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX-1, Z = 1
      index2 = offsetZ + j*offsetY + dimSize-1;
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 1, Z = MAX-1
      index3 = (dimSize-1)*offsetZ + j*offsetY + 1;
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX-1, Z = MAX-1
      index4 = (dimSize-1)*offsetZ + j*offsetY + dimSize-1;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }
    //! X = Y = (0 || MAX)' du/dx = 0, dv/dy = 0
    for (size_t k = 2; k < dimSize-1; k++)
    {
      // X = Y = 0
      uint index1 (k*offsetZ);
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX, Y = 0
      uint index2 (k*offsetZ + dimSize);
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 0, Y = MAX
      uint index3 (k*offsetZ + dimSize*offsetY);
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX, Y = MAX
      uint index4 (k*offsetZ + dimSize*offsetY + dimSize);
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];

      // X = Y = 1
      index1 = k*offsetZ + offsetY + 1;
      triplets0.emplace_back(index1,index1, 1.0);
      // X = MAX-1, Y = 1
      index2 = k*offsetZ + offsetY + dimSize-1;
      triplets0.emplace_back(index2,index2, 1.0);
      // X = 1, Y = MAX-1
      index3 = k*offsetZ + (dimSize-1)*offsetY + 1;
      triplets0.emplace_back(index3,index3, 1.0);
      // X = MAX-1, Y = MAX-1
      index4 = k*offsetZ + (dimSize-1)*offsetY + dimSize-1;
      triplets0.emplace_back(index4,index4, 1.0);
      
      B0[index1] = pExac[index1]; B0[index2] = pExac[index2]; B0[index3] = pExac[index3]; B0[index4] = pExac[index4];
    }
    //-------------------------------------------------------------------------

    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 1
    index = offsetZ + offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = MAX-1, Y = Z = 1
    index = offsetZ + offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = 1, Y = MAX-1
    index = offsetZ + (dimSize-1)*offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = 1, Z = MAX-1
    index = (dimSize-1)*offsetZ + offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = MAX-1, Z = 1
    index = offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = MAX-1, Y = 1
    index = (dimSize-1)*offsetZ + offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // Y = Z = MAX-1, X = 1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + 1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = Z = MAX-1
    index = (dimSize-1)*offsetZ + (dimSize-1)*offsetY + dimSize-1;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];

    // X = Y = Z = 0
    index = 0;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = MAX, Y = Z = 0
    index = dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = 0, Y = MAX
    index = dimSize*offsetY;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = 0, Z = MAX
    index = dimSize*offsetZ;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = MAX, Z = 0
    index = dimSize*offsetY + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Z = MAX, Y = 0
    index = dimSize*offsetZ + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // Y = Z = MAX, X = 0
    index = dimSize*offsetZ + dimSize*offsetY;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    // X = Y = Z = MAX
    index = dimSize*offsetZ + dimSize*offsetY + dimSize;
    triplets0.emplace_back(index,index, 1.0);
    B0[index] = pExac[index];
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // barrier, sync point
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // call Eigen solver for Pressure (idk, does is support concurrency or not)
    //-------------------------------------------------------------------------
    A.setFromTriplets (triplets0.begin(), triplets0.end());
    // biconjugate gradient stabilized algorithm
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double>> solver(A);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Can not build preconditioner" << std::endl;
      std::cout << "Can not build preconditioner" << std::endl;
      return;
    }
    Eigen::VectorXd pHat(vecSize);
    pHat = solver.solveWithGuess(B0, p);
    if (solver.info() != Eigen::Success)
    {
      outputFile << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      std::cout << "Failed to solve the system with Eigen, tick = " << tick << std::endl;
      return;
    }
    //-------------------------------------------------------------------------
    // refresh Pressure values (transfer Eigen-vector to vector)
    //-------------------------------------------------------------------------

    for (size_t i = 0; i < vecSize; ++i) {p0[i] = pHat[i]; p[i] = pHat[i];}

    // residual
    //-------------------------------------------------------------------------
    const double velResidual = velocity_residual(offsetX, offsetY, offsetZ, hX, hY, hZ, tau, dimSize);
    outputResidualFile << std::scientific << velResidual << std::endl;
    //-------------------------------------------------------------------------

    funcOutput(outputFuncFile, "/v1", std::to_string(tick), ".txt", u, params, false);
    funcOutput(outputFuncFile, "/v2", std::to_string(tick), ".txt", v, params, false);
    funcOutput(outputFuncFile, "/v3", std::to_string(tick), ".txt", w, params, false);
    funcOutput(outputFuncFile, "/p", std::to_string(tick), ".txt", p0, params, false);
    tick += 1;
    if ((tick % 100 == 0)) std::cout << "tick: " << tick << '\n';
  }
  for (size_t i = 0; i < vecSize; ++i) p0[i] = p[i];
  std::cout << "final tick: " << tick << '\n';
}
//=================================================================================================

