#include "matrix_io.hpp"
#include "rk4_solver.hpp"
#include <cmath>

//********
//* setup
//********
const std::string dat_dir = "../../dat";
const std::string test_dat_dir = "../../test/dat";
const std::string test_name = "test-rk4_solver-motor";
const std::string dir_prefix = dat_dir + "/" + test_name + "-";
const std::string test_dat_prefix = test_dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string u_arr_fname = "u_arr.dat";
const std::string x_arr_fname = "x_arr.dat";
const std::string x_arr_chk_fname = "x_arr_chk.dat";

constexpr uint_t t_dim = 1e3;
constexpr uint_t x_dim = 3;
constexpr uint_t u_dim = 1;
constexpr real_t t0 = 0;
constexpr real_t x0[x_dim] = {0};
constexpr real_t h = 1. / (t_dim - 1);
// constexpr real_t M_PI = 3.1415926535897932;
#ifdef __USE_SINGLE_PRECISION__
constexpr real_t error_thres = 1e-5;
#else
constexpr real_t error_thres = 1e-6;
#endif
//* Motor equations:
//* dt2__th = -b/J*dt__th + K_t/J*i
//* dt__i = - K_b/L*dt__th - R/L*i + 1/L*e
//*
//* In state-space:
//* x = [th; dt__th; i]
//* u = [e]
//*
//* dt__x = A*x + B*u
//* A = [0, 1,      0;
//*      0, -b/J,   K_t/J;
//*      0, -K_b/L, -R/L];
//* B = [0,
//*      0,
//*      1/L];
//*
//* Crouzet 89830012:
constexpr real_t R = 1.4;      //* [ohm]
constexpr real_t L = 1.7e-3;   //*  [ohm s]
constexpr real_t J = 1.29e-4;  //*  [kg m-2]
constexpr real_t b = 3.92e-4;  //*  [N m s]
constexpr real_t K_t = 6.4e-2; //*  [N m A-1]
constexpr real_t K_b = 6.4e-2; //*  [V s]
constexpr real_t A[x_dim * x_dim] = {0, 1, 0, 0, -b / J, K_t / J, 0, -K_b / L, -R / L};
constexpr real_t B[x_dim * u_dim] = {0, 0, 1 / L};

real_t t = 0;
real_t x[t_dim];
real_t t_arr[t_dim];
real_t x_arr[t_dim * x_dim];
real_t x_arr_chk[t_dim * x_dim];
real_t u_arr[t_dim * u_dim];

//* dt__x = A*x + B*x
void
ode_fun(const real_t, const real_t x[], const uint_t i, real_t OUT_dt__x[])
{
	real_t temp0[x_dim];
	real_t temp1[x_dim];

	const real_t *u = matrix::select_row<u_dim>(i, u_arr);
	matrix::right_multiply<x_dim, x_dim>(A, x, temp0);
	matrix::right_multiply<x_dim, u_dim>(B, u, temp1);
	matrix::sum<x_dim>(temp0, temp1, OUT_dt__x);
}

int
main()
{
	//*****************
	//* read test data
	//*****************
	matrix_io::read<t_dim, u_dim>(test_dat_prefix + u_arr_fname, u_arr);
	matrix_io::read<t_dim, x_dim>(test_dat_prefix + x_arr_chk_fname, x_arr_chk);

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
		const real_t *x = matrix::select_row<x_dim>(i, x_arr);
		const real_t *x_chk = matrix::select_row<x_dim>(i, x_arr_chk);

		for (uint_t j = 0; j < x_dim; ++j) {
			real_t error = std::abs(x[j] - x_chk[j]);
			if (error > max_error) {
				max_error = error;
			}
		}
	}

	//* loop vs cum_loop sanity check
	real_t max_loop_v_cum_loop_error = 0.;
	for (uint_t i = 0; i < x_dim; ++i) {
		real_t error = std::abs(x_arr[x_dim * (t_dim - 1) + i] - x[i]);
		if (error > max_loop_v_cum_loop_error) {
			max_loop_v_cum_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_v_cum_loop_error < error_thres) {
		return 0;
	} else {
		return 1;
	}
}
