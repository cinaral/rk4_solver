#ifndef STEP_HPP_CINARAL_220924_1756
#define STEP_HPP_CINARAL_220924_1756

#include "matrix_op.hpp"
#include "types.hpp"

namespace rk4_solver
{

//* computes the next Runge-Kutta 4th Order step
//*
//* inputs:
//* 1. obj - dynamics object
//* 2. ode_fun - ode function, a member of obj
//* 3. t - time [s]
//* 4. x - [X_DIM] state
//* 5. h - time step [s]
//* 6. i - time index
//*
//* ode_fun can be parametrized using the time (row) index i
template <typename T, uint_t X_DIM>
void
step(T &obj, ode_fun_t<T, X_DIM> ode_fun, const real_t t, const real_t (&x)[X_DIM], const real_t h, const uint_t i,
     real_t (&x_next)[X_DIM])
{
	constexpr real_t rk4_weight_0 = 1. / 6.;
	constexpr real_t rk4_weight_1 = 1. / 3.;
#ifdef __DO_USE_HEAP__ 
	static real_t (*k_0_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (*k_1_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (*k_2_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (*k_3_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (*x_temp_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	real_t (&k_0)[X_DIM] = *k_0_ptr;
	real_t (&k_1)[X_DIM] = *k_1_ptr;
	real_t (&k_2)[X_DIM] = *k_2_ptr;
	real_t (&k_3)[X_DIM] = *k_3_ptr;
	real_t (&x_temp)[X_DIM] = *x_temp_ptr;
	//* "..._ptr"s are of type "real_t(*)[X_DIM]", they point to "real_t[X_DIM]"s which are allocated on the heap,
	//* dereferencing "..._ptr"s gives us rvalue references to "real_t[X_DIM]"s, which can be substituted to "real_t[X_DIM]"s allocated on the stack.
	//* Yes, maybe typedef should have been used more
#else
	static real_t k_0[X_DIM];
	static real_t k_1[X_DIM];
	static real_t k_2[X_DIM];
	static real_t k_3[X_DIM];
	static real_t x_temp[X_DIM];
#endif

	(obj.*ode_fun)(t, x, i, k_0); //* ode_fun(ti, xi)

	//* zero-order hold, i.e. no ODE_FUN(,, i+.5), ODE_FUN(,, i+1,) etc.
	matrix_op::weighted_sum(h / 2, k_0, 1., x, x_temp);
	(obj.*ode_fun)(t + h / 2, x_temp, i, k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

	matrix_op::weighted_sum(h / 2, k_1, 1., x, x_temp);
	(obj.*ode_fun)(t + h / 2, x_temp, i, k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

	matrix_op::weighted_sum(h, k_2, 1., x, x_temp);
	(obj.*ode_fun)(t + h, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

	for (uint_t i = 0; i < X_DIM; ++i) {
		x_next[i] = x[i] +
		    h * (rk4_weight_0 * k_0[i] + rk4_weight_1 * k_1[i] + rk4_weight_1 * k_2[i] + rk4_weight_0 * k_3[i]);
	}
}

} // namespace rk4_solver

#endif