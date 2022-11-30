#ifndef TEST_CONFIG_HPP_CINARAL_221122_1007
#define TEST_CONFIG_HPP_CINARAL_221122_1007

#include "matrix_op.hpp"
#include "matrix_rw.hpp"
#include "rk4_solver.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

using rk4_solver::Real_T;
using rk4_solver::size_t;

namespace test_config
{
//* relative to the directory of the executable
const std::string dat_dir = "../dat";                       //* temporary data directory
const std::string ref_dat_dir = "../../test/reference_dat"; //* reference data directory
const std::string t_arr_fname = "t_arr.dat";                //* file name for the time array
const std::string x_arr_fname = "x_arr.dat";                //* file name for the x array
const std::string x_arr_chk_fname = "x_arr_chk.dat";        //* file name for the reference x array
const std::string u_arr_fname = "u_arr.dat";                //* file name for the u array

template <size_t T_DIM, size_t X_DIM>
Real_T
compute_max_error(const Real_T (&arr)[T_DIM * X_DIM], const Real_T (&arr_chk)[T_DIM * X_DIM])
{
	Real_T max_error = 0.;

	for (size_t i = 0; i < T_DIM; ++i) {
		const Real_T(&a)[X_DIM] = *matrix_op::select_row<T_DIM, X_DIM>(i, arr);
		const Real_T(&a_chk)[X_DIM] = *matrix_op::select_row<T_DIM, X_DIM>(i, arr_chk);

		for (size_t j = 0; j < X_DIM; ++j) {
			const Real_T error = std::abs(a[j] - a_chk[j]);

			if (error > max_error) {
				max_error = error;
			}
		}
	}
	return max_error;
}
} // namespace test_config

#endif