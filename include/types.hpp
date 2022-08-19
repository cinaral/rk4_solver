#ifndef TYPES_HPP_CINARAL_220814_0347
#define TYPES_HPP_CINARAL_220814_0347

using uint_t = unsigned long long int;
#ifdef __USE_SINGLE_PRECISION__
using real_t = float;
#else
using real_t = double;
#endif

using ode_fun_t = void (*)(const real_t t, const real_t x[], const uint_t i, real_t dt__x[]);
using event_fun_t = bool (*)(const real_t t, const real_t x[], const uint_t i, real_t x_plus[]);

#endif