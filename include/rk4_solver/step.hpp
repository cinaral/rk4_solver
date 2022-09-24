#ifndef STEP_HPP_CINARAL_220924_1756
#define STEP_HPP_CINARAL_220924_1756

#include "matrix_op.hpp"
#include "types.hpp"

namespace rk4_solver
{
	
//* computes the next Runge-Kutta 4th Order step
//*
//* inputs:
//* 1. t - time [s]
//* 2. x - [X_DIM] state
//* 3. h - time step [s]
//* 4. i - row index
//*
//* ODE_FUN can be parametrized using i (row_index), but zero-order hold will be used for the parameters during the step
template <uint_t X_DIM, ode_fun_t<X_DIM> ODE_FUN>
void
step(const real_t t, const real_t (&x)[X_DIM], const real_t h, const uint_t i, real_t (&x_next)[X_DIM])
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
	matrix_op::weighted_sum(h / 2, k_0, 1., x, x_temp);
	ODE_FUN(t + h / 2, x_temp, i, k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

	matrix_op::weighted_sum(h / 2, k_1, 1., x, x_temp);
	ODE_FUN(t + h / 2, x_temp, i, k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

	matrix_op::weighted_sum(h, k_2, 1., x, x_temp);
	ODE_FUN(t + h, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

	for (uint_t i = 0; i < X_DIM; ++i) {
		x_next[i] = x[i] +
		    h * (rk4_weight_0 * k_0[i] + rk4_weight_1 * k_1[i] + rk4_weight_1 * k_2[i] + rk4_weight_0 * k_3[i]);
	}
}
} // namespace rk4_solver

#endif