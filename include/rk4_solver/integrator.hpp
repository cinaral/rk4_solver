/*
 * rk4_solver
 *
 * MIT License
 *
 * Copyright (c) 2022 Cinar, A. L.
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

#ifndef INTEGRATOR_HPP_CINARAL_220924_1756
#define INTEGRATOR_HPP_CINARAL_220924_1756

#include "matrix_op.hpp"
#include "types.hpp"

namespace rk4_solver
{
/*
 * Computes the next Runge-Kutta 4th Order step.
 * `ode_fun` can be parametrized using the time (row) index `i`.
 *
 * `step<OPT: X_DIM, T>(obj, ode_fun, t, x, h, (i), OUT:x_next)`
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
template <size_t X_DIM, typename T> class Integrator
{
  public:
	Integrator(T &obj, OdeFun_T<X_DIM, T> ode_fun, const Real_T step_size, const Real_T t0 = 0)
	    : obj(obj), ode_fun(ode_fun), t0(t0), step_size(step_size)
	{
		reset_accumulator();
	}

	void
	step(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&x_next)[X_DIM])
	{
		(obj.*ode_fun)(t, x, i, k_0); //* ode_fun(ti, xi)

		//* zero-order hold, i.e. no ODE_FUN(,, i+.5), ODE_FUN(,, i+1,) etc.
		matrix_op::weighted_sum(step_size / 2, k_0, 1., x, x_temp);
		(obj.*ode_fun)(t + step_size / 2, x_temp, i,
		               k_1); //* ode_fun(ti + h/2, xi + h/2*k_0)

		matrix_op::weighted_sum(step_size / 2, k_1, 1., x, x_temp);
		(obj.*ode_fun)(t + step_size / 2, x_temp, i,
		               k_2); //* ode_fun(ti + h/2, xi + h/2*k_1)

		matrix_op::weighted_sum(step_size, k_2, 1., x, x_temp);
		(obj.*ode_fun)(t + step_size, x_temp, i, k_3); //* ode_fun(ti + h, xi + k_2)

		constexpr Real_T w0 = 1. / 6.;
		constexpr Real_T w1 = 1. / 3.;

		for (size_t i = 0; i < X_DIM; ++i) {
			dx = step_size * (w0 * k_0[i] + w1 * k_1[i] + w1 * k_2[i] + w0 * k_3[i]);
			//* compensated (Kahan) summation, ffast-math might break this
			compensated_dx = dx - accumulator[i];
			x_temp[i] = x[i] + compensated_dx;
			accumulator[i] = (x_temp[i] - x[i]) - compensated_dx;
			x_next[i] = x_temp[i];
		}
	}

	Real_T
	get_step_size() const
	{
		return step_size;
	}

	void
	reset_accumulator()
	{
		for (size_t i = 0; i < X_DIM; ++i) {
			accumulator[i] = 0;
		}
	}

  private:
	T &obj;
	const OdeFun_T<X_DIM, T> ode_fun;
	const Real_T t0;
	const Real_T step_size;
	Real_T dx;
	Real_T compensated_dx;

#ifdef DO_NOT_USE_HEAP
	Real_T k_0[X_DIM];
	Real_T k_1[X_DIM];
	Real_T k_2[X_DIM];
	Real_T k_3[X_DIM];
	Real_T x_temp[X_DIM];
	Real_T accumulator[X_DIM];
#else
	/*
	 * Dereferencing pointers that point to `Real_T[X_DIM]`s which are allocated on the heap, in
	 * order to get rvalue references to the `Real_T[X_DIM]`s.
	 */
	Real_T (&k_0)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T (&k_1)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T (&k_2)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T (&k_3)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T (&x_temp)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
	Real_T (&accumulator)[X_DIM] = *(Real_T(*)[X_DIM]) new Real_T[X_DIM];
#endif
};
} // namespace rk4_solver

#endif
