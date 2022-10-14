#include "matrix_op.hpp"
#include "matrix_rw.hpp"
#include "rk4_solver.hpp"

using Uint_T = rk4_solver::Uint_T;
using Real_T = rk4_solver::Real_T;

//* setup
const std::string dat_dir = "../dat";
const std::string ref_dat_dir = "../../test/reference_dat";
const std::string test_name = "test-rk4_solver-ball";
const std::string dat_prefix = dat_dir + "/" + test_name + "-";
const std::string ref_dat_prefix = ref_dat_dir + "/" + test_name + "-";
const std::string t_arr_fname = "t_arr.dat";
const std::string x_arr_fname = "x_arr.dat";
const std::string x_arr_chk_fname = "x_arr_chk.dat";

constexpr Uint_T sample_freq = 1e4;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0.;
constexpr Real_T t_final = 2.;
constexpr Uint_T t_dim = sample_freq*(t_final - t_init) + 1;
constexpr Uint_T x_dim = 2;
constexpr Real_T x_init[x_dim] = {1., 0.};
constexpr Real_T e_restitution = .75;
constexpr Real_T gravity_const = 9.806;

#ifdef __USE_SINGLE_PRECISION__
constexpr Real_T error_thres = 2e-3;
#else
constexpr Real_T error_thres = 1e-3;
#endif

struct Dynamics {
	//* Ball equations:
	//* dt__x =  [x2; -g]
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const Uint_T, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = -gravity_const;
	}

	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], const Uint_T, Real_T (&x_plus)[x_dim])
	{
		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e_restitution * x[1];
		}
		return false;
	}
};
Dynamics dyn;

int
main()
{
	//* read reference data
	Real_T x_arr_chk[t_dim * x_dim];
	matrix_rw::read<t_dim, x_dim>(ref_dat_prefix + x_arr_chk_fname, x_arr_chk);

	//* test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	rk4_solver::loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t_init, x_init, time_step, &t, x);
	rk4_solver::cum_loop<Dynamics, t_dim, x_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t_init, x_init, time_step, t_arr,
	                                             x_arr);

	//* write test data
	matrix_rw::write<t_dim, 1>(dat_prefix + t_arr_fname, t_arr);
	matrix_rw::write<t_dim, x_dim>(dat_prefix + x_arr_fname, x_arr);

	//* verify
	Real_T max_error = 0.;

	for (Uint_T i = 0; i < t_dim; ++i) {
		const Real_T(&x_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr);
		const Real_T(&x_chk_)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(i, x_arr_chk);

		for (Uint_T j = 0; j < 1; ++j) {
			const Real_T error = std::abs(x_[j] - x_chk_[j]);

			if (error > max_error) {
				max_error = error;
			}
		}
	}

	//* loop vs cum_loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = *matrix_op::select_row<t_dim, x_dim>(t_dim - 1, x_arr);

	for (Uint_T i = 0; i < x_dim; ++i) {
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
