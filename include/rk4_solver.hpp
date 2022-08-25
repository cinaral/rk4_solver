#ifndef RK4_SOLVER_HPP_CINARAL_220814_0005
#define RK4_SOLVER_HPP_CINARAL_220814_0005

//* solves dt__x = ode_fun(t, x(t)), x(0) = x_0 using Runge-Kutta 4th Order Method.

#include "matrix.hpp"

namespace rk4_solver
{

using uint_t = unsigned long long int;
#ifdef __USE_SINGLE_PRECISION__
using real_t = float;
#else
using real_t = double;
#endif

using ode_fun_t = void (*)(const real_t t, const real_t x[], const uint_t i, real_t dt__x[]);
using event_fun_t = bool (*)(const real_t t, const real_t x[], const uint_t i, real_t x_plus[]);

//* computes the next Runge-Kutta 4th Order step
//*
//* inputs:
//* 1. t - time [s]
//* 2. x - [X_DIM] state
//* 3. h - time step [s]
//* 4. i - row index
//*
//* ODE_FUN can be parametrized using i (row_index), but zero-order hold will be used for the parameters during the step
template <ode_fun_t ODE_FUN, uint_t X_DIM>
void
step(const real_t t, const real_t x[], const real_t h, const uint_t i, real_t x_next[])
{
	constexpr real_t rk4_weight_0 = 1. / 6.;
	constexpr real_t rk4_weight_1 = 1. / 3.;

	real_t k_0[X_DIM];
	real_t k_1[X_DIM];
	real_t k_2[X_DIM];
	real_t k_3[X_DIM];
	real_t x_temp[X_DIM];

	ODE_FUN(t, x, i, k_0); //* ode_fun(ti, xi)

	//* zero-order hold, i.e. no ODE_FUN(,, i+.5), ODE_FUN(,, i+1,) etc.
	matrix::weighted_sum<X_DIM>(h / 2, k_0, 1., x, x_temp);
	ODE_FUN(t + h / 2, x_temp, i, k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

	matrix::weighted_sum<X_DIM>(h / 2, k_1, 1., x, x_temp);
	ODE_FUN(t + h / 2, x_temp, i, k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

	matrix::weighted_sum<X_DIM>(h, k_2, 1., x, x_temp);
	ODE_FUN(t + h, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

	for (uint_t i = 0; i < X_DIM; ++i) {
		x_next[i] = x[i] +
		    h * (rk4_weight_0 * k_0[i] + rk4_weight_1 * k_1[i] + rk4_weight_1 * k_2[i] + rk4_weight_0 * k_3[i]);
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void
loop(const real_t t0, const real_t x0[], const real_t h, real_t *t, real_t x[])
{
	matrix::replace_row<X_DIM>(0, x0, x); //* initialize x
	*t = t0;                              //* initialize t

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		step<ODE_FUN, X_DIM>(*t, x, h, i, x); //* update x to the next x

		*t = t0 + (i + 1) * h; //* update t to the next t
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times
//* cumulatively saves all data points into t_arr [T_DIM] and x_arr [T_DIM * X_DIM]
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void
cum_loop(const real_t t0, const real_t x0[], const real_t h, real_t t_arr[], real_t x_arr[])
{
	real_t x[X_DIM] = {0};
	matrix::replace_row<X_DIM>(0, x0, x); //* initialize x
	real_t t = t0;                        //* initialize t

	t_arr[0] = t;
	matrix::replace_row<X_DIM>(0, x, x_arr);

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		step<ODE_FUN, X_DIM>(t, x, h, i, x); //* update x to the next x

		t = t0 + (i + 1) * h; //* update t to the next t

		t_arr[i + 1] = t;
		matrix::replace_row<X_DIM>(i + 1, x, x_arr);
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t
loop(const real_t t0, const real_t x0[], const real_t h, real_t *t, real_t x[])
{
	uint_t i = 0;
	matrix::replace_row<X_DIM>(0, x0, x); //* initialize x
	*t = t0;                              //* initialize t

	//* check for events at the initial condition
	bool stop_flag = EVENT_FUN(*t, x, i, x);

	for (; !stop_flag && i < T_DIM - 1; ++i) {
		step<ODE_FUN, X_DIM>(*t, x, h, i, x); //* update x to the next x

		*t = t0 + (i + 1) * h; //* update t to the next t

		stop_flag = EVENT_FUN(*t, x, i, x);
	}

	return i;
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//* cumulatively saves all data points into t_arr [T_DIM] and x_arr [T_DIM * X_DIM]
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t
cum_loop(const real_t t0, const real_t x0[], const real_t h, real_t t_arr[], real_t x_arr[])
{
	uint_t i = 0;
	real_t x[X_DIM];
	matrix::replace_row<X_DIM>(0, x0, x); //* initialize x
	real_t t = t0;                        //* initialize t

	//* check for events at the initial condition
	bool stop_flag = EVENT_FUN(t, x, i, x);

	t_arr[0] = t;
	matrix::replace_row<X_DIM>(0, x, x_arr);

	for (; !stop_flag && i < T_DIM - 1; ++i) {
		step<ODE_FUN, X_DIM>(t, x, h, i, x); //* update x to the next x

		t = t0 + (i + 1) * h; //* update t to the next t

		stop_flag = EVENT_FUN(t, x, i, x);

		t_arr[i + 1] = t;
		matrix::replace_row<X_DIM>(i + 1, x, x_arr);
	}
	return i;
}

} // namespace rk4_solver

#endif