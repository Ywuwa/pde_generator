#include <fstream>
#include <string>

#include "settings.hpp"
/*! \brief input function
 *  \param[in] filename - the name of input file
 *  \param[in] params   - the link to parameters structure
 *  \return code of operation
 */
int dataInput(const std::string& filename, model_data& params);

/*! \brief grid-function value output
 *  \param[in] path      - path to dir
 *  \param[in] name      - directory name
 *  \param[in] number    - number of file
 *  \param[in] extension - file extension
 *  \param[in] funcValue - function value
 *  \param[in] params    - the link to parameters structure
 *  \param[in] onOffNotification - ON/OFF console text
 *  \return code of operation
 */
void funcOutput(const std::string& path, const std::string& name, const std::string& number, 
	const std::string& extension, const std::vector<double>& funcValue, const model_data& params,
  const bool onOffNotification = false);
  //=================================================================================================