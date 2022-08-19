# ```rk4_solver```: Runge-Kutta 4th Order Solver
[Runge-Kutta 4th Order Method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) ODE Solver with events. This is a header-only C++ library. MATLAB was used to generate reference solutions for testing.

It numerically solves a system of ordinary differential equations (ODE) given as $\dot{\mathbf{x}} = \mathbf{f}(t, \mathbf{x}(t)),\quad \mathbf{x}(0)=\mathbf{x}_0.$

This library has no external dependencies including the standard library. It was written with hard real-time embedded applications in mind, e.g. it does not use dynamic memory allocation, RTTI or exceptions. It can be compiled for ```float```s by enabling the ```__USE_SINGLE_PRECISION__``` compiler flag.

# Installation

Include the headers in ```include/``` into your project. ```matrix_io.hpp``` is optional, it is only used for testing. 

Alternatively, you can use [FetchContent()](https://cmake.org/cmake/help/latest/module/FetchContent.html) in your ```CMakeLists.txt```:
```CMake
FetchContent_Declare(rk4_solver URL https://github.com/cinaral/rk4_solver/releases/download/<RELEASE_TAG>/rk4_solver-src.zip)
FetchContent_MakeAvailable(rk4_solver)
set(rk4_solver_INCLUDE_DIR ${rk4_solver_SOURCE_DIR}/include)
```

Use CTest to test the library before using. Three tests are included:
- Sine wave
- Motor response 
- Bouncing ball
  
Some of them require reference solutions which is provided in this repository without the need for MATLAB, see [Testing](#testing) for details.


# Usage

For a single integration step, call ```rk4_solver::step(...)```:
```Cpp
template <ode_fun_t ODE_FUN, uint_t X_DIM>
void step(const real_t t, const real_t x[], const real_t h, const uint_t i, real_t x_next[]);
```

For an integration loop, call ```rk4_solver::loop(...)```:
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void loop(const real_t t0, const real_t x0[], const real_t h, real_t *t, real_t x[]);
```

If you want to use events with your integration loop, you may do so by using an event function:

**WARNING**: This is a fixed-step size method and therefore the event detection will be approximate within the time step size.
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t loop(const real_t t0, const real_t x0[], const real_t h, real_t *t, real_t x[]);
```

```uint_t```, ```real_t```, ```ode_fun_t``` and ```event_fun_t``` are defined in ```types.hpp```.  

**WARNING:** By default, ```uint_t``` is ```unsigned long long int``` and ```real_t``` is ```double``` (```float``` with ```__USE_SINGLE_PRECISION__``` compiler flag). Their sizes in bytes may be system or compiler dependent. If this is a concern ```types.hpp``` should be modified. 

```ode_fun_t``` and ```event_fun_t``` are function pointers defined as:
```Cpp
using ode_fun_t = void (*)(const real_t t, const real_t x[], const uint_t i, real_t dt__x[]);
using event_fun_t = bool (*)(const real_t t, const real_t x[], const uint_t i, real_t x_plus[]);
```
where ```i < T_DIM``` $\in\mathbb{N}$ is the current time step, i.e. $t[i] = t_0 + i h$ for the time step size $h\in\mathbb{R}$. ```i``` is provided for time index dependent inputs such as pre-computed discrete control input.

You may use cumulative loop functions with or without the event function to save the integration value at every time step:
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void cum_loop(const real_t t0, const real_t x0[], const real_t h, real_t t_arr[], real_t x_arr[]);

template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t cum_loop(const real_t t0, const real_t x0[], const real_t h, real_t t_arr[], real_t x_arr[]);
```


# Examples


## Example 1: Single integration step
```Cpp
#include "rk4_solver.hpp"
#include "types.hpp"
//...
void ode_fun(const real_t, const real_t x[], const uint_t, real_t dt__x[])
{
	dt__x[0] = x[1];
	dt__x[1] = u;
}

int main()
{
	//* integration step
	rk4_solver::step<ode_fun, x_dim>(t, x, h, i, x_next);
	//...
}
```
See [example_step.cpp](./examples/example_step.cpp) for details.


## Example 2: Integration loop
```Cpp
//...
void ode_fun(const real_t t, const real_t[], const uint_t, real_t dt__x[])
{
	dt__x[0] = t;
}

int main()
{
	//* integration loop with cumulatively saved data arrays
	rk4_solver::cum_loop<ode_fun, t_dim, x_dim>(t0, x0, h, t_arr, x_arr);
	//...
}
```
See [example_loop.cpp](./examples/example_loop.cpp) for details.


## Example 3: Events
```Cpp
//...
void ode_fun(const real_t, const real_t x[], const uint_t, real_t dt__x[])
{
	dt__x[0] = x[1];
	dt__x[1] = -g;
}

bool event_fun(const real_t, const real_t x[], const uint_t, real_t x_plus[])
{
	if (x[0] <= 0) {
		x_plus[0] = 0;
		x_plus[1] = -e * x[1];
	}

	//* don't stop the integration
	return false;
}

int main()
{
	//* integration loop with events
	rk4_solver::loop<ode_fun, t_dim, x_dim, event_fun>(t0, x0, h, &t, x);
	//...
}
```
See [example_event.cpp](./examples/example_event.cpp) for details.


# Testing
Reference solutions are required for some tests, which can be found in ```test/dat/```. 

You may need to generate new solutions in order to update the existing tests or to add new tests. [run_all_tests.m](./test/matlab/run_all_tests.m) may be used to generate solutions if you have access to MATLAB. By default, the data files are put in ```dat/```, which you may copy into ```test/dat/``` or use ```scripts/update_test_data.sh```.

The MATLAB tests are optional. The MATLAB scripts under ```test/matlab/``` can also be used to visually inspect the output. 

The ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```.
```MATLAB
writematrix(matrix, file, 'Delimiter', ',');  
matrix = readmatrix(file);  
```
 Each row in file corresponds to a matrix row. For example, if the test reference data is a ```N``` by ```M``` matrix then the ```*.dat``` file should contain ```N``` lines and ```M``` comma-separated numbers in each line.
