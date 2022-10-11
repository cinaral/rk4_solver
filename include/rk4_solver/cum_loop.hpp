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

#ifndef CUM_LOOP_HPP_CINARAL_220924_1755
#define CUM_LOOP_HPP_CINARAL_220924_1755

#include "matrix_op.hpp"
#include "step.hpp"
#include "types.hpp"

namespace rk4_solver
{
//* Loops Runge-Kutta 4th Order step T_DIM times and cumulatively saves all points.
//*
//* cum_loop<T, T_DIM, X_DIM>(obj, ode_fun, t0, x0, h, OUT:t_arr, OUT:x_arr)
//* IN:
//* 1. obj - [T] dynamics object
//* 2. ode_fun - [T::*] ode function, a member of obj
//* 3. t0 - initial time [s]
//* 4. x0 - [X_DIM] initial state
//* 5. h - time step [s]
//* OUT:
//* 6. t_arr - [T_DIM] time array [s]
//* 7. x_arr - [T_DIM * X_DIM] state array
template <typename T, uint_t T_DIM, uint_t X_DIM>
void
cum_loop(T &obj, ode_fun_t<T, X_DIM> ode_fun, const real_t t0, const real_t (&x0)[X_DIM], const real_t h,
         real_t (&t_arr)[T_DIM], real_t (&x_arr)[T_DIM * X_DIM])
{
#ifdef __DO_USE_HEAP__
	static real_t (*x_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (&x)[X_DIM] = *x_ptr;
#else
	real_t x[X_DIM];
#endif
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

//* loops Runge-Kutta 4th Order step T_DIM times or until event_fun returns true and cumulatively saves all points.
//* event_fun can be used to modify x when certain conditions are met.
//*
//* OUT:i = cum_loop<T, T_DIM, X_DIM>(obj, ode_fun, t0, x0, h, OUT:t_arr, OUT:x_arr)
//* IN:
//* 1. obj - [T] dynamics object
//* 2. ode_fun - [T::*] ode function, a member of obj
//* 3. event_fun - [T::*] event function, a member of obj
//* 4. t0 - initial time [s]
//* 5. x0 - [X_DIM] initial state
//* 6. h - time step [s]
//* OUT:
//* 7. i - final index
//* 8. t_arr - [T_DIM] time array [s]
//* 9. x_arr - [T_DIM * X_DIM] state array
template <typename T, uint_t T_DIM, uint_t X_DIM>
uint_t
cum_loop(T &obj, ode_fun_t<T, X_DIM> ode_fun, event_fun_t<T, X_DIM> event_fun, const real_t t0,
         const real_t (&x0)[X_DIM], const real_t h, real_t (&t_arr)[T_DIM], real_t (&x_arr)[T_DIM * X_DIM])
{
	uint_t i = 0;
#ifdef __DO_USE_HEAP__
	static real_t (*x_ptr)[X_DIM] = (real_t(*)[X_DIM])new real_t[X_DIM];
	static real_t (&x)[X_DIM] = *x_ptr;
#else
	real_t x[X_DIM];
#endif
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
