#include "matrix_io.hpp"
#include "rk4_solver.hpp"
#include <cmath>

//********
//* setup
//********
const std::string dat_dir = "../../dat";
const std::string test_name = "test-rk4_solver-sine";
const std::string dir_prefix = dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";

constexpr uint_t t_dim = 1e3;
constexpr uint_t x_dim = 1;
constexpr real_t t0 = 0;
constexpr real_t x0[x_dim] = {0};
constexpr real_t h = 1. / (t_dim - 1);
constexpr real_t M_PI = 3.1415926535897932;
constexpr real_t f = 5.;
#ifdef __USE_SINGLE_PRECISION__
constexpr real_t error_thres = 1e-5;
#else
constexpr real_t error_thres = 1e-9;
#endif

real_t t = 0;
real_t x[x_dim];
real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];
real_t x_arr_chk[t_dim * x_dim];

//* dt__x = f(t, x) = 2*pi*f*cos(t*2*pi*f)
//* x = sin(t*2*pi*f)
void
ode_fun(const real_t t, const real_t[], const uint_t, real_t dt__x[])
{
	dt__x[0] = 2 * M_PI * f * cos(t * 2 * M_PI * f);
}

int
main()
{
	//*******
	//* test
	//*******
	rk4_solver::loop<ode_fun, t_dim, x_dim>(t0, x0, h, &t, x);
	rk4_solver::cum_loop<ode_fun, t_dim, x_dim>(t0, x0, h, t_arr, x_arr);

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

		real_t error = std::abs(x_[0] - sin(t_arr[i] * 2 * M_PI * f));
		if (error > max_error) {
			max_error = error;
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
