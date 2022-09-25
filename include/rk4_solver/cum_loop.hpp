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
//* 1. obj - dynamics object
//* 2. ode_fun - ode function, a member of obj
//* 3. t0 - initial time [s]
//* 4. x0 - [X_DIM] initial state
//* 5. h - time step [s]
template <typename T, uint_t T_DIM, uint_t X_DIM>
void
cum_loop(T &obj, ode_fun_t<T, X_DIM> ode_fun, const real_t t0, const real_t (&x0)[X_DIM], const real_t h,
         real_t (&t_arr)[T_DIM], real_t (&x_arr)[T_DIM * X_DIM])
{
	real_t x[X_DIM] = {0};
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	real_t t = t0;                       //* initialize t

	t_arr[0] = t;
	matrix_op::replace_row<T_DIM>(0, x, x_arr);

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		step<T, X_DIM>(obj, ode_fun, t, x, h, i, x); //* update x to the next x

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
//* 1. obj - dynamics object
//* 2. ode_fun - ode function, a member of obj
//* 3. event_fun - event function, a member of obj
//* 4. t0 - initial time [s]
//* 5. x0 - [X_DIM] initial state
//* 6. h - time step [s]
template <typename T, uint_t T_DIM, uint_t X_DIM>
uint_t
cum_loop(T &obj, ode_fun_t<T, X_DIM> ode_fun, event_fun_t<T, X_DIM> event_fun, const real_t t0,
         const real_t (&x0)[X_DIM], const real_t h, real_t (&t_arr)[T_DIM], real_t (&x_arr)[T_DIM * X_DIM])
{
	uint_t i = 0;
	real_t x[X_DIM];
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	real_t t = t0;                       //* initialize t

	//* check for events at the initial condition
	bool stop_flag = (obj.*event_fun)(t, x, i, x);

	t_arr[0] = t;
	matrix_op::replace_row<T_DIM>(0, x, x_arr);

	for (; !stop_flag && i < T_DIM - 1; ++i) {
		step<T, X_DIM>(obj, ode_fun, t, x, h, i, x); //* update x to the next x

		t = t0 + (i + 1) * h; //* update t to the next t

		stop_flag = (obj.*event_fun)(t, x, i, x);

		t_arr[i + 1] = t;
		matrix_op::replace_row<T_DIM>(i + 1, x, x_arr);
	}
	return i;
}
} // namespace rk4_solver
#endif