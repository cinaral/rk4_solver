#include "test_config.hpp"

//* setup
const std::string test_name = "motor-test";
const std::string dat_prefix = test_config::dat_dir + "/" + test_name + "-";
const std::string ref_dat_prefix = test_config::ref_dat_dir + "/" + test_name + "-";

constexpr size_t sample_freq = 1e3;
constexpr Real_T time_step = 1. / sample_freq;
constexpr Real_T t_init = 0;
constexpr Real_T t_final = 1;
constexpr size_t t_dim = sample_freq * (t_final - t_init) + 1;
constexpr size_t x_dim = 3;
constexpr size_t u_dim = 1;
constexpr Real_T x_init[x_dim] = {0, 0, 0};

constexpr Real_T R = 1.4;      //* [ohm]
constexpr Real_T L = 1.7e-3;   //*  [ohm s]
constexpr Real_T J = 1.29e-4;  //*  [kg m-2]
constexpr Real_T b = 3.92e-4;  //*  [N m s]
constexpr Real_T K_t = 6.4e-2; //*  [N m A-1]
constexpr Real_T K_b = 6.4e-2; //*  [V s]
constexpr Real_T A[x_dim][x_dim] = {{0, 1, 0}, {0, -b / J, K_t / J}, {0, -K_b / L, -R / L}};
constexpr Real_T B[x_dim][u_dim] = {{0}, {0}, {1 / L}};
constexpr Real_T since_ampl = 10; //* input amplitude
constexpr Real_T sine_freq = 10;  //*  input frequency

#ifdef USE_SINGLE_PRECISION
constexpr Real_T error_thres = 1e-5;
#else
constexpr Real_T error_thres = 1e-12;
#endif

struct Dynamics {
	/*
	 * Motor equations:
	 * dt2th = -b/J*dt_th + K_t/J*i
	 * dt_i = - K_b/L*dt_th - R/L*i + 1/L*e
	 *
	 * In state-space:
	 * x = [th; dt_th; i]
	 * u = [e]
	 *
	 * dt_x = A*x + B*u
	 * A = [0, 1,      0;
	 *      0, -b/J,   K_t/J;
	 *      0, -K_b/L, -R/L];
	 * B = [0,
	 *      0,
	 *      1/L];
	 */
	void
	ode_fun(const Real_T t, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		Real_T temp0[x_dim];
		Real_T temp1[x_dim];

		const Real_T u[u_dim] = {
		    static_cast<Real_T>(since_ampl * std::sin(t * 2 * M_PI * sine_freq))};
		matrix_op::right_multiply(A, x, temp0);
		matrix_op::right_multiply(B, u, temp1);
		matrix_op::sum(temp0, temp1, dt_x);
	}
	// Real_T u_arr[t_dim * u_dim];
};
Dynamics dynamics;

int
main()
{
	//* 1. read the reference data
	Real_T x_arr_ref[t_dim][x_dim];
	matrix_rw::read(ref_dat_prefix + test_config::x_arr_ref_fname, x_arr_ref);

	//* 2. test
	Real_T t = 0;
	Real_T x[x_dim];
	Real_T t_arr[1][t_dim];
	Real_T x_arr[t_dim][x_dim];
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);

	rk4_solver::loop<t_dim>(integrator, t_init, x_init, t, x);
	integrator.reset();
	rk4_solver::loop(integrator, t_init, x_init, t_arr[0], x_arr);

	//* 3. write the test data
	matrix_rw::write(dat_prefix + test_config::t_arr_fname, t_arr);
	matrix_rw::write(dat_prefix + test_config::x_arr_fname, x_arr);

	//* 4. verify the results
	Real_T max_error = test_config::compute_max_error(x_arr, x_arr_ref);

	//* loop vs cum_loop sanity check
	Real_T max_loop_error = 0.;
	const Real_T(&x_final)[x_dim] = x_arr[t_dim - 1];

	for (size_t i = 0; i < x_dim; ++i) {
		const Real_T error = std::abs(x_final[i] - x[i]);

		if (error > max_loop_error) {
			max_loop_error = error;
		}
	}

	if (max_error < error_thres && max_loop_error <= std::numeric_limits<Real_T>::epsilon()) {
		return 0;
	} else {
		printf("max_error = %.3g\n", max_error);
		printf("max_loop_error = %.3g\n", max_loop_error);
		return 1;
	}
}
