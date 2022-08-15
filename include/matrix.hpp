#ifndef MATRIX_HPP_CINARAL_220814_0342
#define MATRIX_HPP_CINARAL_220814_0342

//* basic matrix operations

#include "types.hpp"

namespace matrix
{
	
//* selects a row from an N_ROW by M_COL matrix
template <uint_t M_COL>
static void
select_row(const real_t arr[], const uint_t row_idx, real_t sel[])
{
	for (uint_t i = 0; i < M_COL; ++i) {
		sel[i] = arr[row_idx * M_COL + i];
	}
}

//* selects a row from an N_ROW by M_COL matrix
template <uint_t M_COL>
static void
replace_row(const real_t sel[], const uint_t row_idx, real_t arr[])
{
	for (uint_t i = 0; i < M_COL; ++i) {
		arr[row_idx * M_COL + i] = sel[i];
	}
}

//* sums two arrays of DIM size
template <uint_t DIM>
static void
sum(const real_t x[], const real_t y[], real_t OUT_sum[])
{
	for (uint_t i = 0; i < DIM; ++i) {
		OUT_sum[i] = x[i] + y[i];
	}
}

//* scales an array of DIM with a scalar
template <uint_t DIM>
static void
scale(const real_t scale, const real_t x[], real_t OUT_x[])
{
	for (uint_t i = 0; i < DIM; ++i) {
		OUT_x[i] = scale * x[i];
	}
}

//* computes <x, y> for two arrays of DIM size
template <uint_t DIM>
static real_t
dot_product(const real_t x[], const real_t y[])
{
	real_t res = 0;

	for (uint_t i = 0; i < DIM; ++i) {
		res += x[i] * y[i];
	}
	return res;
}

//* computes A*x for an N_ROW by M_COL matrix A
template <uint_t N_ROW, uint_t M_COL>
static void
right_multiply(const real_t A[], const real_t x[], real_t OUT_res[])
{
	real_t row[M_COL];

	for (uint_t i = 0; i < N_ROW; ++i) {
		select_row<M_COL>(A, i, row);
		OUT_res[i] = dot_product<M_COL>(row, x);
	}
}

} // namespace matrix

#endif