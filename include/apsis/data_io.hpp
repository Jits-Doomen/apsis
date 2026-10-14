#pragma once
#include "apsis/orbit.hpp"
#include <string>
#include <vector>
namespace apsis::data {
std::vector<orbit::State> readStateCsv(const std::string& path);
void writeStateCsv(const std::string& path,const std::vector<orbit::State>& states);
void writeXYCsv(const std::string& path,const std::vector<std::pair<double,double>>& data,
                const std::string& x_name="x",const std::string& y_name="y");
}
