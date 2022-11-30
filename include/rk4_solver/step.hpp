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
/*
 * Computes the next Runge-Kutta 4th Order step.
 * `ode_fun` can be parametrized using the time (row) index `i`.
 *
 * `step<OPT: X_DIM, T>(obj, ode_fun, t, x, h, i, OUT:x_next)`
 *
 * 1. `obj`: dynamics object (type `T`)
 * 2. `ode_fun`: ode function, member of `obj` (type `T::*`)
 * 3. `t`: time [s]
 * 4. `x`: state
 * 5. `h`: time step [s]
 * 6. `i`: time index corresponding to `t`
 *
 * OUT:
 *	7. `x_next`: next state
 */
template <size_t X_DIM, typename T>
void
step(T &obj, OdeFun_T<X_DIM, T> ode_fun, const Real_T t, const Real_T (&x)[X_DIM], const Real_T h,
     const size_t i, Real_T (&x_next)[X_DIM])
{
	constexpr Real_T rk4_weight_0 = 1. / 6.;
	constexpr Real_T rk4_weight_1 = 1. / 3.;
#ifdef __DO_USE_HEAP__
	/*
	 * `..._ptr`s are of type `Real_T(*)[X_DIM]`.
	 * They point to `Real_T[X_DIM]`s which are allocated on the heap.
	 * Dereferencing them gives us rvalue references to `Real_T[X_DIM]`s,
	 * which can be substituted for `Real_T[X_DIM]`s allocated on the stack.
	 * (Maybe typedef should be used more.)
	 */
	static Real_T(*k_0_ptr)[X_DIM] = (Real_T(*)[X_DIM]) new Real_T[X_DIM];
	static Real_T(*k_1_ptr)[X_DIM] = (Real_T(*)[X_DIM]) new Real_T[X_DIM];
	static Real_T(*k_2_ptr)[X_DIM] = (Real_T(*)[X_DIM]) new Real_T[X_DIM];
	static Real_T(*k_3_ptr)[X_DIM] = (Real_T(*)[X_DIM]) new Real_T[X_DIM];
	static Real_T(*x_temp_ptr)[X_DIM] = (Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T(&k_0)[X_DIM] = *k_0_ptr;
	Real_T(&k_1)[X_DIM] = *k_1_ptr;
	Real_T(&k_2)[X_DIM] = *k_2_ptr;
	Real_T(&k_3)[X_DIM] = *k_3_ptr;
	Real_T(&x_temp)[X_DIM] = *x_temp_ptr;

#else
	static Real_T k_0[X_DIM];
	static Real_T k_1[X_DIM];
	static Real_T k_2[X_DIM];
	static Real_T k_3[X_DIM];
	static Real_T x_temp[X_DIM];
#endif

	(obj.*ode_fun)(t, x, i, k_0); //* ode_fun(ti, xi)

	//* zero-order hold, i.e. no ODE_FUN(,, i+.5), ODE_FUN(,, i+1,) etc.
	matrix_op::weighted_sum(h / 2, k_0, 1., x, x_temp);
	(obj.*ode_fun)(t + h / 2, x_temp, i, k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

	matrix_op::weighted_sum(h / 2, k_1, 1., x, x_temp);
	(obj.*ode_fun)(t + h / 2, x_temp, i, k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

	matrix_op::weighted_sum(h, k_2, 1., x, x_temp);
	(obj.*ode_fun)(t + h, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

	for (size_t i = 0; i < X_DIM; ++i) {
		x_next[i] = x[i] +
		    h *
			(rk4_weight_0 * k_0[i] + rk4_weight_1 * k_1[i] + rk4_weight_1 * k_2[i] +
		         rk4_weight_0 * k_3[i]);
	}
}
} // namespace rk4_solver

#endif
