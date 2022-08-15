#ifndef MATRIX_io_HPP_CINARAL_220814_17176
#define MATRIX_io_HPP_CINARAL_220814_17176

#include "types.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>


namespace matrix_io
{

const uint_t precision = std::numeric_limits<real_t>::digits10 + 1;
const std::string data_dir = "../../dat";
const std::string delimiter = ",";

template <uint_t N_ROW, uint_t M_COL>
void
write(const std::string file_name, const real_t matrix[])
{
	std::ofstream file;
	file.open(file_name);

	if (file.is_open()) {

		for (uint_t i = 0; i < N_ROW; i++) {

			for (uint_t j = 0; j < M_COL; j++) {
				file << std::setprecision(precision) << std::scientific << matrix[i * M_COL + j];

				if (j < M_COL - 1) {
					file << delimiter;
				}
			}

			if (i < N_ROW) {
				file << std::endl;
			}
		}
	} else {
		std::cout << "Could not open the output file " << file_name << std::endl;
	}
	file.close();
}

template <uint_t N_ROW, uint_t M_COL>
void
read(const std::string file_name, real_t matrix[])
{
	std::ifstream file;
	file.open(file_name);

	if (file.is_open()) {
		std::string line;
		std::string entry;
		// uint_t i = 0;
		// uint_t j = 0;
		size_t str_pos = 0;

		for (uint_t i = 0; i < N_ROW && std::getline(file, line); i++) {
			entry = line.substr(0, line.find(delimiter));

			//* parse the line by splitting at the delimiters
			for (uint_t j = 0; j < M_COL && (str_pos = line.find(delimiter)) != std::string::npos; j++) {
				entry = line.substr(0, str_pos);
				line.erase(0, str_pos + delimiter.length());
				matrix[i * M_COL + j + 1] = std::stof(entry);
			}
			matrix[i * M_COL + M_COL - 1] = std::stof(line);
		}
	} else {
		std::cout << "Could not open the input file " << file_name << std::endl;
	}
	file.close();
}

} // namespace matrix_io

#endif