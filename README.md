
- [1. About rk4_solver](#1-about-rk4solver)
- [2. Installation](#2-installation)
- [3. Usage](#3-usage)
	- [3.1. Single integration step](#31-single-integration-step)
	- [3.2. Integration loop](#32-integration-loop)
	- [3.3. Integration loop with intermediate values](#33-integration-loop-with-intermediate-values)
- [4. Examples](#4-examples)
	- [4.1. Example 1: Single integration step](#41-example-1-single-integration-step)
	- [4.2. Example 2: Integration loop](#42-example-2-integration-loop)
	- [4.3. Example 3: Events](#43-example-3-events)
- [5. Benchmarks](#5-benchmarks)
- [6. Testing](#6-testing)

# 1. About rk4_solver
This is a simple header-only C++ library for solving ordinary differential equations using the [fixed-step 4th order Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods). It offers a rudimentary RK4 solver that is light-weight, fast and hopefully easy to use. The objective is to solve a system of ordinary differential equations (ODEs) given as $\dot{\mathbf{x}} = \mathbf{f}(t, \mathbf{x}(t))$, $\mathbf{x}(0)=\mathbf{x}_0$. 

Some event functionality is included for hybrid dynamics, but it should not be trusted as this is a fixed-step size method.  Kahan summation algorithm is used in the integration to compensate for the accumulation of floating point errors. 

This library and its only dependency [matrix_op](https://github.com/cinaral/matrix_op) was written with embedded systems in mind, and it has been tested on hard real-time embedded systems. It does not use RTTI or exceptions, and you can disable all dynamic heap allocations using the ```DO_NOT_USE_HEAP``` compiler flag, but beware of stack overflows. Use the ```USE_SINGLE_PRECISION``` compiler flag to compile for single precision systems.  

# 2. Installation

Include all or some of the headers in ```include/``` into your project.

Alternatively, you can use [FetchContent()](https://cmake.org/cmake/help/latest/module/FetchContent.html) in your ```CMakeLists.txt```:
```CMake
FetchContent_Declare(rk4_solver URL https://github.com/cinaral/rk4_solver/releases/download/<RELEASE_TAG>/src.zip)
FetchContent_MakeAvailable(rk4_solver)
```

Use CTest to test the library before using. Currently, four tests are available:
1. Integrating a sine function,
2. Solving a first-order system,
3. Solving for the time response of a motor that is driven by a sinusoidal input,
4. Solving for the motion of a bouncing ball (hybrid dynamics).
  
Some of the tests require reference solutions which is included in this repository, see [Testing](#testing) for details. For the minimum examples, see ```examples/```.

# 3. Usage

You must first create an integrator object:
```Cpp
rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step, t_init = 0);
```

ODE function should be a member of ```Dynamics```, which is a class you will provide for the dynamics of your system. The ODE function must be of type ```OdeFun_T```:
```Cpp
//* dx/dt = ode_fun(t, x)
void Dynamics::ode_fun(t, x, OUT: dt_x); 
```
If you need an event function, you can use an event function of type ```EventFun_T```. It is fine to use the same object for both the integrator and the event function:
```Cpp
//* x_plus = event_fun(t, x)
void Events::check(t, x, OUT: x_plus); 
```
**WARNING**: This is a fixed-step size method and therefore the event detection will be approximate within the time step size.
```OdeFun_T``` and ```EventFun_T``` are function pointers defined in [types.hpp](include/rk4_solver/types.hpp):  
```Cpp
using OdeFun_T = void (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&dt_x)[X_DIM]);
using EventFun_T = bool (T::*)(const Real_T t, const Real_T (&x)[X_DIM], const size_t i, Real_T (&x_plus)[X_DIM]);
```
**WARNING:** By default, ```Real_T``` is ```double```. Use ```USE_SINGLE_PRECISION``` compiler flag to set to ```float```.

## 3.1. Single integration step
Call ```Integrator::step(...)```:
```Cpp
void
step(
	const Real_T t, 
	const Real_T (&x)[X_DIM], 
	Real_T (&x_next)[X_DIM]
);
```

## 3.2. Integration loop
Call ```rk4_solver::loop(...)``` as follows if you want to discard all intermediate values except the final state:
```Cpp
void
loop(
	Integrator integrator, 
	OPTIONAL: Event event,
	const Real_T t_init,
	const Real_T (&x0)[X_DIM], 
	Real_T &t, 
	Real_T (&x)[X_DIM]
);
```
## 3.3. Integration loop with intermediate values
Otherwise, you can use the overloaded version of ```loop(...)``` that takes an array for time and another array for state history to save all intermediate values of the solution:
```Cpp
void
loop(
	Integrator integrator, 
	OPTIONAL: Event event,
	const Real_T t_init,
	const Real_T (&x0)[X_DIM], 
	Real_T (&t)[T_DIM], 
	Real_T (&x)[T_DIM * X_DIM]
);
```
Where ```i < T_DIM``` is the current time step. 

# 4. Examples

## 4.1. Example 1: Single integration step
```Cpp
#include "rk4_solver/integrator.hpp"
//...
	Dynamics::ode_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&dt_x)[x_dim])
	{
		dt_x[0] = x[1];
		dt_x[1] = u;
	}
//...
Dynamics dynamics;

int main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	Real_T t;
	Real_T x[x_dim];
	integrator.step(t_init, x_init, t, x); //* integration step
	//...
}
```
See [step-example.cpp](./examples/step-example.cpp) for details.


## 4.2. Example 2: Integration loop
```Cpp
#include "rk4_solver/cum_loop.hpp"
//...
int main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	Real_T t_arr[t_dim];
	Real_T x_arr[t_dim * x_dim];
	//* save the intermediate values
	rk4_solver::loop<t_dim>(integrator, t_init, x_init, t_arr, x_arr);
	//...
}
```
See [loop-example.cpp](./examples/loop-example.cpp) for details.


## 4.3. Example 3: Events
```Cpp
#include "rk4_solver.hpp"
//...
struct Dynamics {
	//...
	bool
	event_fun(const Real_T, const Real_T (&x)[x_dim], Real_T (&x_plus)[x_dim])
	{
		bool did_occur = false;

		if (x[0] <= 0) {
			x_plus[0] = 0;
			x_plus[1] = -e_restitution * x[1];
			did_occur = true;
		}
		return did_occur;
	}
};
Dynamics dynamics;

int
main()
{
	rk4_solver::Integrator<x_dim, Dynamics> integrator(dynamics, &Dynamics::ode_fun, time_step);
	rk4_solver::Event<x_dim, Dynamics> event(dynamics, &Dynamics::event_fun);
	Real_T t;
	Real_T x[x_dim];
	//* discard the intermediate values and stop at first event
	rk4_solver::loop<t_dim>(integrator, event, t_init, x_init, t, x, true);
	//...
}
```
See [example-event.cpp](./examples/event-example.cpp) for details.

# 5. Benchmarks

There are two benchmark tests: 
1. A step integration loop without final time, and intermediate values are discarded.
2. A cumulative integration loop with final time, and intermediate values are saved.

The benchmark test is a 3rd order linear system compiled using g++ with ```-O3``` optimization level. A desktop Intel i7-9700K at 3.60 GHz processor with 32 GB of memory was used to obtain the following sample benchmarks: 

|                                                  Flags | Loop (million steps per second) | Cumulative Loop (million steps per second) |
| -----------------------------------------------------: | :-----------------------------: | :----------------------------------------: |
|                                       None *(Default)* |              18.5               |                    28.6                    |
|                             ```USE_SINGLE_PRECISION``` |              18.5               |                    30.5                    |
|                                  ```DO_NOT_USE_HEAP``` |              18.6               |                    28.2                    |
| ```DO_NOT_USE_HEAP``` *and* ```USE_SINGLE_PRECISION``` |              18.6               |                    30.7                    |


Using the ```USE_SINGLE_PRECISION``` or the ```DO_NOT_USE_HEAP``` flag does not affect the performance meaningfully.

**WARNING**: Your stack can easily overflow for large problems with the ```DO_NOT_USE_HEAP``` flag and should be avoided if it is not necessary.

# 6. Testing
Reference solutions are required for some tests, which are included in ```test/dat/```. [matrix_rw](https://github.com/cinaral/matrix_rw) library is used to read and write from file, the ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```.

You may need to generate new reference solutions to update the existing tests or add new tests. [run_test.m](./test/matlab/run_test.m) can be used to (re)generate the reference solutions if you have MATLAB. The MATLAB tests in ```test/matlab/``` are optional, but they can be used to visually verify the results. 

