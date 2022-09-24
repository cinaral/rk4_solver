#ifndef CUM_LOOP_HPP_CINARAL_220924_1755
#define CUM_LOOP_HPP_CINARAL_220924_1755

#include "matrix_op.hpp"
#include "step.hpp"
#include "types.hpp"

namespace rk4_solver
{
//* loops Runge-Kutta 4th Order step T_DIM times
//* cumulatively saves all data points into t_arr [T_DIM] and x_arr [T_DIM * X_DIM]
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <uint_t T_DIM, uint_t X_DIM, ode_fun_t<X_DIM> ODE_FUN>
void
cum_loop(const real_t t0, const real_t (&x0)[X_DIM], const real_t h, real_t (&t_arr)[T_DIM],
         real_t (&x_arr)[T_DIM * X_DIM])
{
	real_t x[X_DIM] = {0};
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	real_t t = t0;                       //* initialize t

	t_arr[0] = t;
	matrix_op::replace_row<T_DIM>(0, x, x_arr);

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		step<X_DIM, ODE_FUN>(t, x, h, i, x); //* update x to the next x

		t = t0 + (i + 1) * h; //* update t to the next t

		t_arr[i + 1] = t;
		matrix_op::replace_row<T_DIM>(i + 1, x, x_arr);
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//* cumulatively saves all data points into t_arr [T_DIM] and x_arr [T_DIM * X_DIM]
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <uint_t T_DIM, uint_t X_DIM, ode_fun_t<X_DIM> ODE_FUN, event_fun_t<X_DIM> EVENT_FUN>
uint_t
cum_loop(const real_t t0, const real_t (&x0)[X_DIM], const real_t h, real_t (&t_arr)[T_DIM],
         real_t (&x_arr)[T_DIM * X_DIM])
{
	uint_t i = 0;
	real_t x[X_DIM];
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	real_t t = t0;                       //* initialize t

	//* check for events at the initial condition
	bool stop_flag = EVENT_FUN(t, x, i, x);

	t_arr[0] = t;
	matrix_op::replace_row<T_DIM>(0, x, x_arr);

	for (; !stop_flag && i < T_DIM - 1; ++i) {
		step<X_DIM, ODE_FUN>(t, x, h, i, x); //* update x to the next x

		t = t0 + (i + 1) * h; //* update t to the next t

		stop_flag = EVENT_FUN(t, x, i, x);

		t_arr[i + 1] = t;
		matrix_op::replace_row<T_DIM>(i + 1, x, x_arr);
	}
	return i;
}
} // namespace rk4_solver
#endif