#include "matrix_op.hpp"
#include "matrix_rw.hpp"
#include "rk4_solver.hpp"

using size_t = rk4_solver::size_t;
using Real_T = rk4_solver::Real_T;

//* setup
const std::string dat_dir = "../dat";
const std::string ref_dat_dir = "../../test/reference_dat";
const std::string test_name = "test-rk4_solver-motor";
const std::string dat_prefix = dat_dir + "/" + test_name + "-";
const std::string ref_dat_prefix = ref_dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string u_arr_fname = "u_arr.dat";
const std::string x_arr_fname = "x_arr.dat";
const std::string x_arr_chk_fname = "x_arr_chk.dat";

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0;
constexpr Real_T t_final = 1;
constexpr size_t t_dim = sample_freq*(t_final - t_init) + 1;
constexpr size_t x_dim = 3;
constexpr size_t u_dim = 1;
constexpr Real_T x_init[x_dim] = {0,0,0};

constexpr Real_T R = 1.4;      //* [ohm]
constexpr Real_T L = 1.7e-3;   //*  [ohm s]
constexpr Real_T J = 1.29e-4;  //*  [kg m-2]
constexpr Real_T b = 3.92e-4;  //*  [N m s]
constexpr Real_T K_t = 6.4e-2; //*  [N m A-1]
constexpr Real_T K_b = 6.4e-2; //*  [V s]
constexpr Real_T A[x_dim * x_dim] = {0, 1, 0, 0, -b / J, K_t / J, 0, -K_b / L, -R / L};
constexpr Real_T B[x_dim * u_dim] = {0, 0, 1 / L};

#ifdef __USE_SINGLE_PRECISION__
constexpr Real_T error_thres = 1e-5;
#else
constexpr Real_T error_thres = 1e-13;
#endif

struct Dynamics {
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
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t i, Real_T (&dt__x)[x_dim])
	{
		Real_T temp0[x_dim];
		Real_T temp1[x_dim];

		const Real_T(&u)[u_dim] = *matrix_op::select_row<t_dim, u_dim>(i, u_arr);
		matrix_op::right_multiply<x_dim, x_dim>(A, x, temp0);
		matrix_op::right_multiply<x_dim, u_dim>(B, u, temp1);
		matrix_op::sum<x_dim>(temp0, temp1, dt__x);
	}
	Real_T u_arr[t_dim * u_dim];
};
Dynamics dyn;

int
main()
{
	//* read reference data
	Real_T x_arr_chk[t_dim * x_dim];
	matrix_rw::read<t_dim, u_dim>(ref_dat_prefix + u_arr_fname, dyn.u_arr);
	matrix_rw::read<t_dim, x_dim>(ref_dat_prefix + x_arr_chk_fname, x_arr_chk);

	//* test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, &t, x);
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, t_init, x_init, time_step, t_arr, x_arr);

	//* write test data
	matrix_rw::write<t_dim, 1>(dat_prefix + t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + x_arr_fname, x_arr);

	//* verify
	Real_T max_error = 0.;

	for (size_t i = 0; i < t_dim; ++i) {
		const Real_T(&x_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr);
		const Real_T(&x_chk_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr_chk);

		for (size_t j = 0; j < x_dim; ++j) {
			const Real_T error = std::abs(x_[j] - x_chk_[j]);

			if (error > max_error) {
				max_error = error;
			}
		}
	}

	//* loop vs cum_loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(t_dim - 1, x_arr);

	for (size_t i = 0; i < x_dim; ++i) {
		const Real_T error = std::abs(x_final[i] - x[i]);

		if (error > max_loop_error) {
			max_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_error < std::numeric_limits<Real_T>::epsilon()) {
		return 0;
	} else {
		return 1;
	}
}
