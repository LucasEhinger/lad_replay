#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void change_gem_param(const std::string& filename, double x0, double y0, double z0,
                      double x1, double y1, double z1, double theta) {
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    return;
  }

  std::ostringstream updatedContent;
  std::string line;

  while (std::getline(infile, line)) {
    if (line.find("lgem_m0_position") == 0) {
      updatedContent << "lgem_m0_position = " << x0 << ", " << y0 << ", " << z0 << std::endl;
    } else if (line.find("lgem_m1_position") == 0) {
      updatedContent << "lgem_m1_position = " << x1 << ", " << y1 << ", " << z1 << std::endl;
    } else if (line.find("lgem_angle") == 0) {
      updatedContent << "lgem_angle = " << theta << std::endl;
    }
    else {
      updatedContent << line << std::endl;
    }
  }

  infile.close();

  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not write to file " << filename << std::endl;
    return;
  }

  outfile << updatedContent.str();
  outfile.close();
}

int main() {
  const std::string filename = "lgem.param";
  double x = 1.5, y = 0.0, z = 95.571;

  // change_gem_param(filename, x, y, z);

  return 0;
}