#include "matrix_io.hpp"
#include "rk4_solver.hpp"
#include "types.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>

constexpr real_t M_PI = 3.1415926535897932;

//********
//* setup
//********
const std::string dat_dir = "../../dat";
const std::string test_dat_dir = "../../test/dat";
const std::string test_name = "test-rk4_solver-sine";
const std::string dir_prefix = dat_dir + "/" + test_name + "-";
const std::string test_dat_prefix = test_dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";
const std::string x_arr_chk_fname = "x_arr_chk.dat";

const uint_t t_dim = 1e3;
const uint_t x_dim = 1;
const real_t error_thres = 1e-6;

real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];
real_t x_arr_chk[t_dim * x_dim];

//* dt__x = f(t, x) = 2*pi*cos(t)
void
ode_fun(const real_t t, const real_t[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = 2 * M_PI * cos(t * 2 * M_PI);
}

int
main()
{
	//*****************
	//* read test data
	//*****************
	matrix_io::read<t_dim, 1>(test_dat_prefix + t_arr_fname, t_arr);
	matrix_io::read<t_dim, 1>(test_dat_prefix + x_arr_chk_fname, x_arr_chk);

	//*******
	//* test
	//*******
	rk4_solver::loop<ode_fun, t_dim, x_dim>(t_arr, x_arr);

	//******************
	//* write test data
	//******************
	matrix_io::write<t_dim, x_dim>(dir_prefix + x_arr_fname, x_arr);

	//*********
	//* verify
	//*********
	real_t max_error = 0.;

	for (uint_t i = 0; i < t_dim; i++) {

		real_t x[x_dim];
		real_t x_chk[x_dim];

		matrix::select_row<x_dim>(x_arr, i, x);
		matrix::select_row<x_dim>(x_arr_chk, i, x_chk);

		for (uint_t j = 0; j < x_dim; j++) {
			real_t error = std::abs(x[j] - x_chk[j]);
			if (error > max_error) {
				max_error = error;
			}
		}
	}

	if (max_error < error_thres) {
		return 0;
	} else {
		return 1;
	}
}
