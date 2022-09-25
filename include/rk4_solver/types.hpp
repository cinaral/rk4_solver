#ifndef TYPES_HPP_CINARAL_220924_1756
#define TYPES_HPP_CINARAL_220924_1756

namespace rk4_solver
{
using uint_t = unsigned long long int;
#ifdef __USE_SINGLE_PRECISION__
using real_t = float;
#else
using real_t = double;
#endif

template <typename T, uint_t X_DIM>
using ode_fun_t = void (T::*)(const real_t t, const real_t (&x)[X_DIM], const uint_t i, real_t (&dt__x)[X_DIM]);

template <typename T, uint_t X_DIM>
using event_fun_t = bool (T::*)(const real_t t, const real_t (&x)[X_DIM], const uint_t i, real_t (&x_plus)[X_DIM]);

} // namespace rk4_solver

#endif