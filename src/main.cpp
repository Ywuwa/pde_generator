#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/compute_flow.hpp"
#include <chrono>

//=================================================================================================
int main(int argc, char* argv[]) {
  model_data model;
  std::filesystem::path full_path(std::filesystem::absolute(argv[0]));
  std::filesystem::path dir = full_path.parent_path().parent_path();
  model.PATH = dir.string();
  std::cout << model.PATH << '\n';
  const std::string inputFile = model.PATH + model.PATH_config;
  int code = dataInput(inputFile, model);
  if (code < 0) return 0;
  model.show();

  // function init
  const uint vecSize ( (model.domainPartition + 1) * 
                        (model.domainPartition + 1) * 
                        (model.domainPartition + 1) );
  std::vector<double> velX(vecSize, 0.0);
  initialConditions(velX, 0, model);
  std::vector<double> velY(vecSize, 0.0);
  initialConditions(velY, 1, model);
  std::vector<double> velZ(vecSize, 0.0);
  initialConditions(velZ, 2, model);
  std::vector<double> pressure(vecSize, 0.0);
  initialConditions(pressure, 3, model);

  const std::string outputFuncFile = model.PATH;
  funcOutput(outputFuncFile, "/v1", std::to_string(0), ".txt", velX, model, true);
  funcOutput(outputFuncFile, "/v2", std::to_string(0), ".txt", velY, model, true);
  funcOutput(outputFuncFile, "/v3", std::to_string(0), ".txt", velZ, model, true);
  funcOutput(outputFuncFile, "/p", std::to_string(0), ".txt", pressure, model, true);

  auto start = std::chrono::high_resolution_clock::now();
  // flow computation
  compute_cube(model, velX, velY, velZ, pressure);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Время выполнения: " << elapsed.count() << " секунд" << std::endl;
  
  std::ofstream outputFile(model.PATH + model.PATH_log, std::ios::app);
  if (outputFile.is_open()) outputFile << "execution is finished" << '\n';
}
//=================================================================================================
