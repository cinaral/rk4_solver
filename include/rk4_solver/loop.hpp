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

#ifndef LOOP_HPP_CINARAL_220924_1755
#define LOOP_HPP_CINARAL_220924_1755

#include "event.hpp"
#include "integrator.hpp"
#include "matrix_op.hpp"
#include "types.hpp"

namespace rk4_solver
{
/*
 * Loops Runge-Kutta 4th Order step `T_DIM` times.
 *
 * `loop<T_DIM, OPT: X_DIM, T>(obj, ode_fun, t0, x0, h, OUT:t, OUT:x)`
 *
 * 1. `obj`: dynamics object (type `T`)
 * 2. `ode_fun`: ode function, member of `obj` (type `T::*`)
 * 3. `t0`: initial time [s]
 * 4. `x0`: initial state
 * 5. `h`: time step [s]
 *
 * OUT:
 * 6. `t`: final time [s]
 * 7. `x`: final state
 */
template <size_t T_DIM, size_t X_DIM, typename T>
void
loop(Integrator<X_DIM, T> integrator, const Real_T t0, const Real_T (&x0)[X_DIM], Real_T *t,
     Real_T (&x)[X_DIM])
{
	static const Real_T h = integrator.get_step_size();

	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	*t = t0;                             //* initialize t

	for (size_t i = 0; i < T_DIM - 1; ++i) {
		integrator.step(*t, x, i, x); //* update x to the next x
		*t = t0 + (i + 1) * h;        //* update t to the next t
	}
}

/*
 * Loops Runge-Kutta 4th Order step `T_DIM` times or until event_fun returns true.
 * `event_fun` can be used to modify x when certain conditions are met.
 *
 * `i = loop<T_DIM, OPT:  X_DIM, T>(obj, ode_fun, event_fun, t0, x0, h, OUT:t, OUT:x)`
 * `i`: final index
 *
 * 1. `obj`: dynamics object (type `T`)
 * 2. `ode_fun`: ode function, member of `obj` (type `T::*`)
 * 3. `event_fun`: event function, member of `obj` (type `T::*`)
 * 4. `t0`: initial time [s]
 * 5. `x0`: initial state
 * 6. `h`: time step [s]
 *
 * OUT:
 * 6. `t`: final/event time [s]
 * 7. `x`: final/event state
 */
template <size_t T_DIM, size_t X_DIM, typename T>
size_t
loop(Integrator<X_DIM, T> integrator, Event<X_DIM, T> event, const Real_T t0,
     const Real_T (&x0)[X_DIM], Real_T *t, Real_T (&x)[X_DIM])
{
	static const Real_T h = integrator.get_step_size();

	size_t i = 0;
	matrix_op::replace_row<1>(0, x0, x); //* initialize x
	*t = t0;                             //* initialize t

	for (; i < T_DIM - 1 && !event.check(*t, x, i, x); ++i) {
		integrator.step(*t, x, i, x); //* update x to the next x
		*t = t0 + (i + 1) * h;        //* update t to the next t
	}
	return i;
}
} // namespace rk4_solver
#endif
