# ```rk4_solver```: Runge-Kutta 4th Order Solver
Runge-Kutta 4th Order Method ODE Solver with events. This is a header-only library.

It numerically solves a system of ordinary differential equations (ODE) given as $\dot{\mathbf{x}} = \mathbf{f}(t, \mathbf{x}(t)),\quad \mathbf{x}(0)=\mathbf{x}_0.$


# Usage

Include the headers in ```include/``` into your project. ```matrix_io.hpp``` is optional.

For a single integration step, call ```rk4_solver::step(...)```:
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void step(const real_t t, const real_t x[], const uint_t i, const real_t h, real_t OUT_x_next[])
```

For an integration loop, call ```rk4_solver::loop(...)```:
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM>
void loop(const real_t t_arr[], real_t x_arr[])
```

If you want to use events with your integration loop, you may do so by using an event function:
```Cpp
template <ode_fun_t ODE_FUN, uint_t T_DIM, uint_t X_DIM, event_fun_t EVENT_FUN>
uint_t loop(const real_t t_arr[], real_t x_arr[])
```

```uint_t```, ```real_t```, ```ode_fun_t``` and ```event_fun_t``` are defined in ```types.hpp```.  By default, ```uint_t``` is a 32-bit unsigned integer and ```real_t``` is double. ```ode_fun_t``` and ```event_fun_t``` are function pointers defined as:
```Cpp
using ode_fun_t = void (*)(const real_t t, const real_t x[], const uint_t i, real_t OUT_dt__x[]);
using event_fun_t = bool (*)(const uint_t i, real_t x[]);
```
# Examples

## Example 1: Single integration step
```Cpp
#include "rk4_solver.hpp"
#include "types.hpp"
//...
void ode_fun(const real_t, const real_t x[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = x[1];
	OUT_dt__x[1] = u;
}

int main()
{
	rk4_solver::step<ode_fun, x_dim>(0, x, 0, h, OUT_x_next);
	//...
}
```
See [example_step.cpp](./examples/example_step.cpp) for details.


## Example 2: Integration loop
```Cpp
//...
void ode_fun(const real_t t, const real_t[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = t;
}

int main()
{
	rk4_solver::loop<ode_fun, t_dim, x_dim>(t_arr, x_arr);
	//...
}
```
See [example_loop.cpp](./examples/example_loop.cpp) for details.

## Example 3: Events
<!--See ```examples/example_loop.cpp``` for details.-->
```Cpp
//...
void ode_fun(const real_t t, const real_t[], const uint_t, real_t OUT_dt__x[])
{
	OUT_dt__x[0] = t;
}

int main()
{
	rk4_solver::loop<ode_fun, t_dim, x_dim, event_fun>(t_arr, x_arr);
	//...
}
```


# To do

1. Bouncing ball test
2. Change array access to pointer iteration: (for int* p = ARR; p < DIM; ++p) *p = x;
3. Functions to operate on forward iterators, no more temp arrays, same idea as above but "p = BEGIN; p < END"
4. Reduce logical structures, function calls

# Testing
Reference data is required for some of the tests, which can be found in ```test/dat/```. 

You may need to generate new data in order to update the existing tests or to add new tests. [run_all_tests.m](./test/matlab/run_all_tests.m) can be used to generate reference data if you have access to MATLAB. By default the data files are put in ```dat/```, which you may copy into ```test/dat/```. 

The ```*.dat``` files are comma and newline delimited. If you have access to MATLAB, the formatting is compatible with ```writematrix``` and ```readmatrix```.
```MATLAB
writematrix(matrix, file, 'Delimiter', ',');  
matrix = readmatrix(file);  
```
 Each row in file corresponds to a matrix row. For example, if the test reference data is an ```N``` by ```M``` matrix then the ```*.dat``` file should contain ```N``` rows and ```M``` numbers at each row.
