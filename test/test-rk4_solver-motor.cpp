//#include "matrix_io.hpp"
//#include "rk4_solver.hpp"
//#include "types.hpp"

//#include <fstream>
//#include <iomanip>
//#include <limits>
////#include <cmath>

//const std::string data_dir = "../../dat";
//const std::string test_name = "test";
//const std::string t_arr_fname = "t_arr.dat";
//const std::string x_arr_fname = "x_arr.dat";
//const std::string u_arr_fname = "u_arr.dat";
//const std::string dir_prefix = data_dir + "/" + test_name + "-";

//const int t_dim = 101;

//constexpr uint_t x_dim = 3;
//constexpr uint_t u_dim = 1;
//constexpr real_t R = 1.4; //* [ohm]
//constexpr real_t L = 1.7e-3; //*  [ohm s]
//constexpr real_t J = 1.29e-4; //*  [kg m-2]
//constexpr real_t b = 3.92e-4; //*  [N m s]
//constexpr real_t K_t = 6.4e-2; //*  [N m A-1]
//constexpr real_t K_b = 6.4e-2; //*  [V s]
//constexpr real_t A[x_dim * x_dim] = {0, 1, 0, 0, -b / J, K_t / J, 0, -K_b / L, -R / L};
//constexpr real_t B[x_dim * u_dim] = {0, 0, 1 / L};

//float t;
//float sine;

//real_t u_arr[t_dim * u_dim];

//void
//ode_fun(const real_t, const real_t x[], const uint_t i, real_t OUT_dt__x[])
//{
//	//* dt__x = A*x + B*x
//	real_t temp0[x_dim];
//	real_t temp1[x_dim];
//	real_t u[u_dim];
//	matrix::select_row<x_dim>(u_arr, i, u);
//	matrix::right_multiply<x_dim, x_dim>(A, x, temp0);
//	matrix::right_multiply<x_dim, u_dim>(B, u, temp1);
//	matrix::sum<x_dim>(temp0, temp1, OUT_dt__x);
//}

int
main()
{
	return 0;
}
