#ifndef RK4_SOLVER_HPP_CINARAL_220814_0005
#define RK4_SOLVER_HPP_CINARAL_220814_0005

//* solves dt__x = ode_fun(t, x(t)), x(0) = x_0 using Runge-Kutta 4th Order Method.

#include "matrix.hpp"
#include "types.hpp"

//#ifdef __AVX__
//  #include <immintrin.h>
//#else
//  #warning AVX is not available. Code will not compile!
//#endif

namespace rk4_solver
{

//* computes the next Runge-Kutta 4th Order step
//*
//* inputs:
//* 1. t - time
//* 2. x[X_DIM] - state
//* 3. i - row index
//*
//* ODE_FUN can be parametrized using i (row_index), but zero-order hold will be used for the parameters during the step
template <ode_fun_t ODE_FUN, uint_t X_DIM>
void
step(const real_t t, const real_t x[], const uint_t i, const real_t h, real_t OUT_x_next[])
{
	constexpr real_t rk4_weight_0 = 1. / 6.;
	constexpr real_t rk4_weight_1 = 1. / 3.;
	const real_t t_next_0 = t + h / 2;
	const real_t t_next_1 = t + h;

	real_t k_0[X_DIM];
	real_t k_1[X_DIM];
	real_t k_2[X_DIM];
	real_t k_3[X_DIM];
	real_t x_temp[X_DIM];

	ODE_FUN(t, x, i, k_0); //* ode_fun(ti, xi)

	//* zero-order hold, i.e. no ODE_FUN(,, i+.5), ODE_FUN(,, i+1,) etc.
	matrix::weighted_sum<X_DIM>(h / 2, k_0, 1., x, x_temp);
	ODE_FUN(t_next_0, x_temp, i, k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

	matrix::weighted_sum<X_DIM>(h / 2, k_1, 1., x, x_temp);
	ODE_FUN(t_next_0, x_temp, i, k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

	matrix::weighted_sum<X_DIM>(h, k_2, 1., x, x_temp);
	ODE_FUN(t_next_1, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

	for (uint_t i = 0; i < X_DIM; ++i) {
		OUT_x_next[i] = x[i] +
		    h * (rk4_weight_0 * k_0[i] + rk4_weight_1 * k_1[i] + rk4_weight_1 * k_2[i] + rk4_weight_0 * k_3[i]);
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times
//*
//* inputs:
//* 1. h - time step
//* 2. t - time
//* 3. x - state
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void
loop(const real_t h, real_t *t, real_t x[])
{
	real_t x_next[X_DIM] = {0};

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		*t = i * h;

		step<ODE_FUN, X_DIM>(*t, x, i, h, x_next);
		matrix::replace_row<X_DIM>(x_next, 0, x);
	}

	*t = (T_DIM - 1) * h;
}

//* loops Runge-Kutta 4th Order step T_DIM times
//* cumulatively saves all data points
//*
//* inputs:
//* 1. t_arr[T_DIM] - time array
//* 2. x_arr[T_DIM x X_DIM] - state array
//*
//* the initial condition is supplied in the first row of x_arr:
//* t_arr[0] = t_0;
//*	x_arr[0, X_DIM] = x_0
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void
cum_loop(const real_t t_arr[], real_t x_arr[])
{
	real_t x[X_DIM] = {0};
	real_t x_next[X_DIM] = {0};

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		const real_t t = t_arr[i];
		const real_t h = t_arr[i + 1] - t_arr[i];

		matrix::select_row<X_DIM>(x_arr, i, x);
		step<ODE_FUN, X_DIM>(t, x, i, h, x_next);
		matrix::replace_row<X_DIM>(x_next, i + 1, x_arr);
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//*
//* inputs:
//* 1. h - time step
//* 2. t - time
//* 3. x - state
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t
loop(const real_t h, real_t *t, real_t x[])
{
	uint_t i = 0;
	bool stop_flag = EVENT_FUN(i, x);

	for (; i < T_DIM - 1 && !stop_flag; ++i) {
		*t = i * h;

		step<ODE_FUN, X_DIM>(*t, x, i, h, x);
		stop_flag = EVENT_FUN(i, x);
	}

	*t = (T_DIM - 1) * h;

	return i;
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//* cumulatively saves all data points
//*
//* inputs:
//* 1. t_arr[T_DIM] - time array
//* 2. x_arr[T_DIM x X_DIM] - state array
//*
//* outputs:
//* 1. i - event row index
//*
//* the initial condition is supplied in the first row of x_arr:
//* t_arr[0] = t_0;
//*	x_arr[0, X_DIM] = x_0
//*
//*
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t
cum_loop(const real_t t_arr[], real_t x_arr[])
{
	uint_t i = 0;
	real_t x[X_DIM];

	//* check for events at the initial condition
	matrix::select_row<X_DIM>(x_arr, 0, x);
	bool stop_flag = EVENT_FUN(i, x);
	matrix::replace_row<X_DIM>(0, x, x_arr);

	for (; i < T_DIM - 1 && !stop_flag; ++i) {
		const real_t t = t_arr[i];
		const real_t h = t_arr[i + 1] - i;

		real_t x_next[X_DIM];

		matrix::select_row<X_DIM>(x_arr, i, x);
		step<ODE_FUN, T_DIM, X_DIM>(t, x, i, h, x_next);
		stop_flag = EVENT_FUN(i, x_next);
		matrix::replace_row<X_DIM>(x_next, i + 1, x_arr);
	}

	return i;
}

} // namespace rk4_solver

#endif