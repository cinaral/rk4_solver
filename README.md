# ```rk4_solver```: Runge-Kutta 4th Order Solver
Simple [Runge-Kutta 4th Order Method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) ODE Solver with events. This is a header-only C++ library. MATLAB was used to generate reference solutions for testing.

It numerically solves a system of ordinary differential equations (ODE) given as $\dot{\mathbf{x}} = \mathbf{f}(t, \mathbf{x}(t)),\quad \mathbf{x}(0)=\mathbf{x}_0$.

This library and its only dependency [matrix_op](https://github.com/cinaral/matrix_op) was written with hard real-time embedded applications in mind (excluding [Testing](#testing)), e.g. it does not use dynamic memory allocation, RTTI or exceptions. It can be compiled for ```float```s by enabling the ```__USE_SINGLE_PRECISION__``` compiler flag.


# Installation

Include all or some of the headers in ```include/``` into your project.

Alternatively, you can use [FetchContent()](https://cmake.org/cmake/help/latest/module/FetchContent.html) in your ```CMakeLists.txt```:
```CMake
FetchContent_Declare(rk4_solver URL https://github.com/cinaral/rk4_solver/releases/download/<RELEASE_TAG>/src.zip)
FetchContent_MakeAvailable(rk4_solver)
```

Use CTest to test the library before using. There are three included tests.
- Sine wave
- Motor response 
- Bouncing ball
  
Some of them require reference solutions which is provided in this repository without the need for MATLAB, see [Testing](#testing) for details.


# Usage

For a single integration step, call ```rk4_solver::step(...)```:
```Cpp
template <size_t T_DIM, size_t X_DIM, typename T>
void 
loop(
	T &obj,
	OdeFun_T<T, X_DIM> ode_fun,
	const Real_T t0,
	const Real_T (&x0)[X_DIM],
	const Real_T h,
	Real_T *t,
	Real_T (&x)[X_DIM]
);
```

For an integration loop, call ```rk4_solver::loop(...)```:
```Cpp
template <size_t T_DIM, size_t X_DIM, typename T>
void
loop(
	T &obj,
	OdeFun_T<T, X_DIM> ode_fun,
	const Real_T t0,
	const Real_T (&x0)[X_DIM],
	const Real_T h,
	Real_T *t,
	Real_T (&x)[X_DIM]
);
```

If you want to use events with your integration loop, you may do so by using an event function:

**WARNING**: This is a fixed-step size method and therefore the event detection will be approximate within the time step size.
```Cpp
template <size_t T_DIM, size_t X_DIM, typename T>
size_t
loop(
	T &obj,
	OdeFun_T<T, X_DIM> ode_fun,
	EventFun_T<T, X_DIM> event_fun,
	const Real_T t0,
	const Real_T (&x0)[X_DIM],
	const Real_T h,
	Real_T *t,
	Real_T (&x)[X_DIM]
);
```

**WARNING:** By default, ```Real_T``` is ```double```. Use ```__USE_SINGLE_PRECISION__``` compiler flag to set to ```float```.

```OdeFun_T``` and ```EventFun_T``` are function pointers defined in [types.hpp](include/rk4_solver/types.hpp):  
```Cpp
template <size_t X_DIM, typename T>
using OdeFun_T = void (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&dt__x)[X_DIM]);

template <size_t X_DIM, typename T>
using EventFun_T = bool (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&x_plus)[X_DIM]);

```
Where ```i < T_DIM``` is the current time step for the time step size ```h```. ```i``` is provided for time index dependent inputs such as pre-computed discrete control input.

You may use cumulative loop functions with or without the event function to save the integration value at every time step:
```Cpp
template <size_t T_DIM, size_t X_DIM, typename T>
void
cum_loop(	
	T &obj, 	
	OdeFun_T<T, X_DIM> ode_fun, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	const Real_T h, 
	Real_T (&t_arr)[T_DIM], 
	Real_T (&x_arr)[T_DIM * X_DIM]
);

template <size_t T_DIM, size_t X_DIM, typename T>
size_t
cum_loop(	
	T &obj, 
	OdeFun_T<T, X_DIM> ode_fun, 
	EventFun_T<T, X_DIM> event_fun, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	const Real_T h, 
	Real_T (&t_arr)[T_DIM], 
	Real_T (&x_arr)[T_DIM * X_DIM]
);
```


# Examples

## Example 1: Single integration step
```Cpp
#include "rk4_solver/step.hpp"
//...
void Dynamics::ode_fun(const Real_T, const Real_T x[], const size_t, Real_T dt__x[])
{
	dt__x[0] = x[1];
	dt__x[1] = u;
}
//...
Dynamics dyn;

int main()
{
	//* integration step
	rk4_solver::step(dyn, &Dynamics::ode_fun, t, x, h, i, x_next);
	//...
}
```
See [step-example.cpp](./examples/step-example.cpp) for details.


## Example 2: Integration loop
```Cpp
#include "rk4_solver/cum_loop.hpp"
//...
void Dynamics::ode_fun(const Real_T t, const Real_T[], const size_t, Real_T dt__x[])
{
	dt__x[0] = t;
}
//...
Dynamics dyn;

int main()
{
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<t_dim>(dyn, &Dynamics::ode_fun, t0, x0, h, t_arr, x_arr);
	//...
}
```
See [loop-example.cpp](./examples/loop-example.cpp) for details.


## Example 3: Events
```Cpp
#include "rk4_solver.hpp"
//...
struct Dynamics {
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt__x)[x_dim])
	{
		dt__x[0] = x[1];
		dt__x[1] = -g;
	}
	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&x_plus)[x_dim])
	{
		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e * x[1];
		}
		//* don't stop the integration
		return false;
	}
};
Dynamics dyn;

int
main()
{
	//* integration loop with events
	rk4_solver::loop<t_dim>(dyn, &Dynamics::ode_fun, &Dynamics::event_fun, t0, x0, h, &t, x);
	//...
}
```
See [example-event.cpp](./examples/event-example.cpp) for details.


# Testing
Reference solutions are required for some tests, which can be found in ```test/dat/```. 

[matrix_rw](https://github.com/cinaral/matrix_rw) library is used for testing, the ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```. 

You may need to generate new solutions in order to update the existing tests or to add new tests. [run_test.m](./test/matlab/run_test.m) may be used to generate solutions if you have access to MATLAB. By default, the data files are put in ```dat/```, which you may copy into ```test/dat/``` or use ```scripts/update_test_data.sh```.

The MATLAB tests are optional. The MATLAB scripts under ```test/matlab/``` can also be used to visually inspect the output. 
