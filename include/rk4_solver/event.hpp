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

#ifndef EVENT_HPP_CINARAL_230321_1039
#define EVENT_HPP_CINARAL_230321_1039

#include "types.hpp"

namespace rk4_solver
{
template <size_t X_DIM, typename T> class Event
{
  public:
	Event(T &obj, EventFun_T<X_DIM, T> event_fun) : obj(obj), event_fun(event_fun)
	{
	}

	int
	check(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&x_next)[X_DIM])
	{
		return (obj.*event_fun)(t, x, i, x_next);
	}

  private:
	T &obj;
	EventFun_T<X_DIM, T> event_fun;
};
} // namespace rk4_solver
#endif