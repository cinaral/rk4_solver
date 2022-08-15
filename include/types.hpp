#ifndef TYPES_HPP_CINARAL_220814_0347
#define TYPES_HPP_CINARAL_220814_0347

#include <cstdint>

using uint_t = uint32_t;
#ifdef __USE_SINGLE_PRECISION__
using real_t = float;
#else
using real_t = double;
#endif

using ode_fun_t = void (*)(const real_t t, const real_t x[], const uint_t i, real_t OUT_dt__x[]);
using event_fun_t = bool (*)(const uint_t i, real_t x[]);

#endif