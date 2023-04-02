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
 * 1. `integrator`: integrator object
 * 2. `t_init`: initial time [s]
 * 3. `x_init`: initial state
 *
 * OUT:
 * 4. `t`: final time [s]
 * 5. `x`: final state
 */
template <size_t T_DIM, size_t X_DIM, typename T>
void
loop(Integrator<X_DIM, T> integrator, const Real_T &t_init, const Real_T (&x_init)[X_DIM],
     Real_T &t, Real_T (&x)[X_DIM])
{
	t = t_init; //* initialize t

	for (size_t i = 0; i < X_DIM; ++i) {
		x[i] = x_init[i]; //* initialize x
	}

	for (size_t i = 0; i < T_DIM - 1; ++i) {
		integrator.step(t, x, t, x); //* update t, x to the next t, x
	}
}

/*
 * Loops Runge-Kutta 4th Order step `T_DIM` times and cumulatively saves the results
 *
 * 1. `integrator`: integrator object
 * 2. `t_init`: initial time [s]
 * 3. `x_init`: initial state
 *
 * OUT:
 * 4. `t_arr`: time history
 * 5. `x_arr`: state history
 */
template <size_t T_DIM, size_t X_DIM, typename T>
void
loop(Integrator<X_DIM, T> integrator, const Real_T &t_init, const Real_T (&x_init)[X_DIM],
     Real_T (&t_arr)[T_DIM], Real_T (&x_arr)[T_DIM][X_DIM])
{
	t_arr[0] = t_init;                               //* initialize t
	matrix_op::replace_row<T_DIM>(0, x_init, x_arr); //* initialize x

	Real_T t_next;
	Real_T x_next[X_DIM];

	for (size_t i = 0; i < T_DIM - 1; ++i) {
		const Real_T t = t_arr[i];
		const Real_T(&x)[X_DIM] = x_arr[i];
		integrator.step(t, x, t_next, x_next); //* update t, x to the next t, x
		t_arr[i + 1] = t_next;
		matrix_op::replace_row<T_DIM>(i + 1, x_next, x_arr);
	}
}

/*
 * Loops Runge-Kutta 4th Order step `T_DIM` times or until event_fun returns true.
 * `event_fun` can be used to modify x when certain conditions are met.
 *
 * 1. `integrator`: integrator object
 * 2. `integrator`: event object
 * 3. `t_init`: initial time [s]
 * 4. `x_init`: initial state
 *
 * OUT:
 * 5. `t`: final time
 * 6. `x`: final state
 */
template <size_t T_DIM, size_t X_DIM, typename T>
size_t
loop(Integrator<X_DIM, T> integrator, Event<X_DIM, T> event, const Real_T &t_init,
     const Real_T (&x_init)[X_DIM], Real_T &t, Real_T (&x)[X_DIM], bool halt_on_event = false)
{
	t = t_init; //* initialize t

	for (size_t i = 0; i < X_DIM; ++i) {
		x[i] = x_init[i]; //* initialize x
	}
	Real_T x_plus[X_DIM];

	for (size_t i = 0; i < T_DIM - 1; ++i) {
		if (event.check(t, x, x_plus)) {
			for (size_t j = 0; j < X_DIM; ++j) {
				x[j] = x_plus[j];
			}

			if (halt_on_event) {
				break;
			}
		}
		integrator.step(t, x, t, x); //* update t, x to the next t, x
	}
	return integrator.get_step_count();
}

/*
 * Loops Runge-Kutta 4th Order step `T_DIM` times or until event_fun returns true and cumulatively
 * saves all points. `event_fun` can be used to modify x when certain conditions are met.
 *
 * 1. `integrator`: integrator object
 * 2. `integrator`: event object
 * 3. `t_init`: initial time [s]
 * 4. `x_init`: initial state
 *
 * OUT:
 * 5. `t_arr`: time history
 * 6. `x_arr`: state history
 */
template <size_t T_DIM, size_t X_DIM, typename T>
size_t
loop(Integrator<X_DIM, T> integrator, Event<X_DIM, T> event, const Real_T &t_init,
     const Real_T (&x_init)[X_DIM], Real_T (&t_arr)[T_DIM], Real_T (&x_arr)[T_DIM][X_DIM],
     bool halt_on_event = false)
{
	t_arr[0] = t_init;                               //* initialize t
	matrix_op::replace_row<T_DIM>(0, x_init, x_arr); //* initialize x
	Real_T t_next;
	Real_T x_next[X_DIM];
	Real_T x_plus[X_DIM];

	for (size_t i = 0; i < T_DIM - 1; ++i) {
		const Real_T t = t_arr[i];
		const Real_T(&x)[X_DIM] = x_arr[i];

		if (event.check(t, x, x_plus)) {
			integrator.step(t, x_plus, t_next, x_next);
			t_arr[i + 1] = t_next;
			matrix_op::replace_row<T_DIM>(i + 1, x_next, x_arr);

			if (halt_on_event) {
				break;
			}
		} else {
			integrator.step(t, x, t_next, x_next); //* update t, x to the next t, x
			t_arr[i + 1] = t_next;
			matrix_op::replace_row<T_DIM>(i + 1, x_next, x_arr);
		}
	}
	return integrator.get_step_count();
}

} // namespace rk4_solver
#endif
