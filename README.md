
- [1. About ```rk4_solver```](#1-about-rk4solver)
- [2. Installation](#2-installation)
- [3. Usage](#3-usage)
- [4. Testing](#4-testing)
- [5. Examples](#5-examples)
	- [5.1. Example 1: Single integration step](#51-example-1-single-integration-step)
	- [5.2. Example 2: Integration loop](#52-example-2-integration-loop)
	- [5.3. Example 3: Events](#53-example-3-events)
- [6. Benchmarks](#6-benchmarks)
	- [6.1. Discussion](#61-discussion)

# 1. About ```rk4_solver```
This is a simple header-only C++ library for solving ordinary differential equations (ODEs) with events using the [Runge-Kutta 4th Order Method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods). The objective is to solve a system of ODEs given as $\dot{\mathbf{x}} = \mathbf{f}(t, \mathbf{x}(t))$, $\mathbf{x}(0)=\mathbf{x}_0$. 

This library and its only dependency [matrix_op](https://github.com/cinaral/matrix_op) was written with hard real-time embedded applications in mind (excluding [Testing](#testing)), e.g. it does not use dynamic memory allocation, RTTI or exceptions. It can be compiled for ```float```s by enabling the ```USE_SINGLE_PRECISION``` compiler flag.

# 2. Installation

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


# 3. Usage

For integration, create an integrator object first:
```Cpp
rk4_solver::Integrator<X_DIM, T> integrator(T &obj, OdeFun_T<X_DIM, T> ode_fun, const Real_T step_size, const Real_T t0 = 0);
```

For a single integration step, call ```Integrator::step(...)```:
```Cpp
void
step(
	const Real_T t, 
	const Real_T (&x)[X_DIM], 
	const size_t i, 
	Real_T (&x_next)[X_DIM]
);
```
If you want to loop, call ```rk4_solver::loop(...)```:
```Cpp
void
loop(
	Integrator integrator, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	Real_T *t, 
	Real_T (&x)[X_DIM]
);
```

If you need events, you want an event object:
```Cpp
	rk4_solver::Event<X_DIM, T> event(T &obj, EventFun_T<X_DIM, T> event_fun);
```
**WARNING**: This is a fixed-step size method and therefore the event detection will be approximate within the time step size.
```Cpp
size_t
loop(
	Integrator integrator, 
	Event event, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	Real_T *t, 
	Real_T (&x)[X_DIM]
);
```

**WARNING:** By default, ```Real_T``` is ```double```. Use ```USE_SINGLE_PRECISION``` compiler flag to set to ```float```.

```OdeFun_T``` and ```EventFun_T``` are function pointers defined in [types.hpp](include/rk4_solver/types.hpp):  
```Cpp
using OdeFun_T = void (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&dt_x)[X_DIM]);
using EventFun_T = bool (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&x_plus)[X_DIM]);
```
Where ```i < T_DIM``` is the current time step for the time step size ```h```. ```i``` is provided for time index dependent inputs such as pre-computed discrete control input.

You may use cumulative integration loop functions with or without the event function to save the integration value at every time step:
```Cpp
void
cum_loop(
	Integrator integrator, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	Real_T (&t_arr)[T_DIM], 
	Real_T (&x_arr)[T_DIM * X_DIM]
);

size_t
cum_loop(
	Integrator integrator, 
	Event event, 
	const Real_T t0,
	const Real_T (&x0)[X_DIM], 
	Real_T (&t_arr)[T_DIM], 
	Real_T (&x_arr)[T_DIM * X_DIM]
);
```

# 4. Testing
Reference solutions are required for some tests, which can be found in ```test/dat/```. 

[matrix_rw](https://github.com/cinaral/matrix_rw) library is used for testing, the ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```. 

You may need to generate new solutions in order to update the existing tests or to add new tests. [run_test.m](./test/matlab/run_test.m) may be used to generate solutions if you have access to MATLAB. By default, the data files are put in ```dat/```, which you may copy into ```test/dat/``` or use ```scripts/update_test_data.sh```.

The MATLAB tests are optional, and the scripts under ```test/matlab/``` can be used to visualize the results. 

# 5. Examples

## 5.1. Example 1: Single integration step
```Cpp
#include "rk4_solver/integrator.hpp"
//...
void Dynamics::ode_fun(const Real_T, const Real_T x[], const size_t, Real_T dt_x[])
{
	dt_x[0] = x[1];
	dt_x[1] = u;
}
//...
Dynamics dynamics;

return 0;
int main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, h);
	Real_T x_next[x_dim];
	integrator.step(t, x, i, x_next); //* integration step
	//...
}
```
See [step-example.cpp](./examples/step-example.cpp) for details.


## 5.2. Example 2: Integration loop
```Cpp
#include "rk4_solver/cum_loop.hpp"
//...
void Dynamics::ode_fun(const Real_T t, const Real_T[], const size_t, Real_T dt_x[])
{
	dt_x[0] = t;
}
//...
Dynamics dynamics;

int main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, h);
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<t_dim>(integrator, t0, x0, t_arr, x_arr);
	//...
}
```
See [cum_loop-example.cpp](./examples/cum_loop-example.cpp) for details.


## 5.3. Example 3: Events
```Cpp
#include "rk4_solver.hpp"
//...
struct Dynamics {
	void
	ode_fun(const Real_T, const Real_T (&x)[x_dim], const size_t, Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = -g;
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
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, h);
	rk4_solver::Event<x_dim, Dynamics> event(dynamics, &Dynamics::event_fun);
	Real_T t;
	Real_T x[x_dim];
	//* integration loop with events
	rk4_solver::loop<t_dim>(integrator, event, t0, x0, &t, x);
	//...
}
```
See [example-event.cpp](./examples/event-example.cpp) for details.

# 6. Benchmarks

There are two benchmark tests: 
1. A step integration loop without final time, and intermediate values are discarded.
2. A cumulative integration loop with final time, and intermediate values are saved.

The benchmark test is a 3rd order linear system compiled using g++ with ```-O3``` optimization level. An Intel i7-9700K at 3.60 GHz processor with 32 GB of memory was used to obtain the following results: 

|                                                  Flags | Loop (million steps per second) | Cumulative Loop (million steps per second) |
| -----------------------------------------------------: | :-----------------------------: | :----------------------------------------: |
|                                       None *(Default)* |              18.6               |                    33.3                    |
|                             ```USE_SINGLE_PRECISION``` |              18.6               |                    36.1                    |
|                                  ```DO_NOT_USE_HEAP``` |              18.6               |                    25.6                    |
| ```DO_NOT_USE_HEAP``` *and* ```USE_SINGLE_PRECISION``` |              22.3               |                    28.6                    |



## 6.1. Discussion
1. Using the ```USE_SINGLE_PRECISION``` flag to use single-precision floats does not affect the performance.
2. Using the ```DO_NOT_USE_HEAP``` flag seems to benefit loops without final time, but it hurts cumulative loops. However, what is happening here could be simply due to different memory layouts, so take these results with a grain of salt.
3. If the problem size can fit the stack size, then using the ```DO_NOT_USE_HEAP``` flag to disable heap allocation can provide a performance boost to loops without final time.
**WARNING**: Your stack can easily overflow for large problems with the ```DO_NOT_USE_HEAP``` flag.
