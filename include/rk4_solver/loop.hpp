#ifndef LOOP_HPP_CINARAL_220924_1755
#define LOOP_HPP_CINARAL_220924_1755

#include "matrix_op.hpp"
#include "step.hpp"
#include "types.hpp"

namespace rk4_solver
{
//* loops Runge-Kutta 4th Order step T_DIM times
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <uint_t T_DIM, uint_t X_DIM, ode_fun_t<X_DIM> ODE_FUN>
void
loop(const real_t t0, const real_t (&x0)[X_DIM], const real_t h, real_t *t, real_t (&x)[X_DIM])
{
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	*t = t0;                             //* initialize t

	for (uint_t i = 0; i < T_DIM - 1; ++i) {
		step<X_DIM, ODE_FUN>(*t, x, h, i, x); //* update x to the next x

		*t = t0 + (i + 1) * h; //* update t to the next t
	}
}

//* loops Runge-Kutta 4th Order step T_DIM times or until EVENT_FUN returns true.
//* EVENT_FUN can be used to modify x when certain conditions are met.
//*
//* inputs:
//* 1. t0 - initial time [s]
//* 2. x0 - [X_DIM] initial state
//* 3. h - time step [s]
template <uint_t T_DIM, uint_t X_DIM, ode_fun_t<X_DIM> ODE_FUN, event_fun_t<X_DIM> EVENT_FUN>
uint_t
loop(const real_t t0, const real_t (&x0)[X_DIM], const real_t h, real_t *t, real_t (&x)[X_DIM])
{
	uint_t i = 0;
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	*t = t0;                             //* initialize t

	//* check for events at the initial condition
	bool stop_flag = EVENT_FUN(*t, x, i, x);

	for (; !stop_flag && i < T_DIM - 1; ++i) {
		step<X_DIM, ODE_FUN>(*t, x, h, i, x); //* update x to the next x

		*t = t0 + (i + 1) * h; //* update t to the next t

		stop_flag = EVENT_FUN(*t, x, i, x);
	}

	return i;
}
} // namespace rk4_solver
#endif