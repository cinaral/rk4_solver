//* requires input data for verification
//* test__rk4_solver__ball.m can generate it in ./dat if you have MATLAB (see README.md)
//* then copy to ./test/dat or use ./scripts/update_test_data.sh

#include "matrix_io.hpp"
#include "rk4_solver.hpp"

//********
//* setup
//********
const std::string dat_dir = "../../dat";
const std::string test_dat_dir = "../../test/dat";
const std::string test_name = "test-rk4_solver-ball";
const std::string dir_prefix = dat_dir + "/" + test_name + "-";
const std::string test_dat_prefix = test_dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";
const std::string x_arr_chk_fname = "x_arr_chk.dat";

constexpr uint_t t_dim = 1e5;
constexpr uint_t x_dim = 2;
constexpr real_t t0 = 0.;
constexpr real_t tf = 2.;
constexpr real_t x0[x_dim] = {1., 0.};
constexpr real_t h = tf / (t_dim - 1);
#ifdef __USE_SINGLE_PRECISION__
constexpr real_t error_thres = 1e-3;
#else
constexpr real_t error_thres = 1e-4;
#endif

//* Ball equations:
//* dt__x =  [x2; -g]
constexpr real_t e = .75;
constexpr real_t g = 9.806;

real_t t = 0;
real_t x[x_dim];
real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];
real_t x_arr_chk[t_dim * x_dim];

void
ode_fun(const real_t, const real_t x[], const uint_t, real_t dt__x[])
{
	dt__x[0] = x[1];
	dt__x[1] = -g;
}
bool
event_fun(const real_t, const real_t x[], const uint_t, real_t x_plus[])
{
	if (x[0] <= 0) {
		x_plus[0] = 0;
		x_plus[1] = -e * x[1];
	}
	return false;
}

int
main()
{
	//*****************
	//* read test data
	//*****************
	matrix_io::read<t_dim, x_dim>(test_dat_prefix + x_arr_chk_fname, x_arr_chk);

	//*******
	//* test
	//*******
	rk4_solver::loop<ode_fun, t_dim, x_dim, event_fun>(t0, x0, h, &t, x);
	rk4_solver::cum_loop<ode_fun, t_dim, x_dim, event_fun>(t0, x0, h, t_arr, x_arr);

	//******************
	//* write test data
	//******************
	matrix_io::write<t_dim, 1>(dir_prefix + t_arr_fname, t_arr);
	matrix_io::write<t_dim, x_dim>(dir_prefix + x_arr_fname, x_arr);

	//*********
	//* verify
	//*********
	real_t max_error = 0.;

	for (uint_t i = 0; i < t_dim; ++i) {
		const real_t *x_ = matrix::select_row<x_dim>(i, x_arr);
		const real_t *x_chk_ = matrix::select_row<x_dim>(i, x_arr_chk);

		for (uint_t j = 0; j < 1; ++j) {
			real_t error = std::abs(x_[j] - x_chk_[j]);
			if (error > max_error) {
				max_error = error;
			}
		}
	}

	//* loop vs cum_loop sanity check
	real_t max_loop_vs_cum_error = 0.;
	for (uint_t i = 0; i < x_dim; ++i) {
		real_t error = std::abs(x_arr[x_dim * (t_dim - 1) + i] - x[i]);
		if (error > max_loop_vs_cum_error) {
			max_loop_vs_cum_error = error;
		}
	}

	if (max_error < error_thres && max_loop_vs_cum_error < error_thres) {
		return 0;
	} else {
		return 1;
	}
}
