/*
 * rk4_solver
 *  
 * MIT License
 * 
 * Copyright (c) 2022 cinaral
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef STEP_HPP_CINARAL_220924_1756
#define STEP_HPP_CINARAL_220924_1756

#include "matrix_op.hpp"
#include "types.hpp"

namespace rk4_solver
{

//* Computes the next Runge-Kutta 4th Order step.
//* ode_fun can be parametrized using the time (row) index i.
//*
//* step<T, X_DIM>(obj, ode_fun, t, x, h, i, OUT:x_next)
//* IN:
//* 1. obj - [T] dynamics object
//* 2. ode_fun - [T::*] ode function, a member of obj
//* 3. t - time [s]
//* 4. x - [X_DIM] state
//* 5. h - time step [s]
//* 6. i - time index
//* OUT:
//*	7. x_next - [X_DIM] next state
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
