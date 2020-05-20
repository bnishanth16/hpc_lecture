// Name: Nishanth Baskaran
// Student ID: 19M15017
// HPSC Assignment-L4
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main()
{
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for (int i = 0; i < N; i++)
  {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }

  float j_id[N];
  __m256 xv = _mm256_load_ps(x);
  __m256 yv = _mm256_load_ps(y);
  __m256 mv = _mm256_load_ps(m);

  for (int i = 0; i < N; i++)
  {
    __m256 iv = _mm256_set1_ps(i);
    __m256 xiv = _mm256_set1_ps(x[i]);
    __m256 yiv = _mm256_set1_ps(y[i]);

    for (int j = 0; j < N; j++) j_id[j] = j;
      
    __m256 jv = _mm256_load_ps(j_id);
    __m256 mask = _mm256_cmp_ps(jv, iv, _CMP_NEQ_OQ);

    __m256 xjv = _mm256_setzero_ps();
    __m256 yjv = _mm256_setzero_ps();
    __m256 mjv = _mm256_setzero_ps();

    xjv = _mm256_blendv_ps(xjv, xv, mask);
    yjv = _mm256_blendv_ps(yjv, yv, mask);
    mjv = _mm256_blendv_ps(mjv, mv, mask);

    // rx = x[i] - x[j];
    __m256 rxv = _mm256_sub_ps(xiv, xjv);
    // ry = y[i] - y[j];
    __m256 ryv = _mm256_sub_ps(yiv, yjv);

    // rx * rx
    __m256 rxv_square = _mm256_mul_ps(rxv, rxv);
    // ry * ry
    __m256 ryv_square = _mm256_mul_ps(ryv, ryv);

    // (rx * rx + ry * ry)
    __m256 rv_add = _mm256_add_ps(rxv_square, ryv_square);
    // (1/r) = 1/sqrt(rx * rx + ry * ry);
    __m256 rv = _mm256_rsqrt_ps(rv_add); //reciprocal of square root

    // 1/ (r * r * r)
    __m256 rv_cube = _mm256_mul_ps(rv,_mm256_mul_ps(rv, rv));

    // fx[i] -= rx * m[j] / (r * r * r);
    __m256 fxvi = -_mm256_mul_ps(_mm256_mul_ps(rxv, mjv), rv_cube);
    __m256 fxv = _mm256_permute2f128_ps(fxvi, fxv, 1);
    fxv = _mm256_add_ps(fxv, fxvi);
    fxv = _mm256_hadd_ps(_mm256_hadd_ps(fxv, fxv), _mm256_hadd_ps(fxv, fxv));

    // fy[i] -= ry * m[j] / (r * r * r);
    __m256 fyvi = -_mm256_mul_ps(_mm256_mul_ps(ryv, mjv), rv_cube);
    __m256 fyv = _mm256_permute2f128_ps(fyvi, fyv, 1);
    fyv = _mm256_add_ps(fyv, fyvi);
    fyv = _mm256_hadd_ps(_mm256_hadd_ps(fyv, fyv), _mm256_hadd_ps(fyv, fyv));

    // storing values back to vectors
    _mm256_store_ps(fx, fxv);
    _mm256_store_ps(fy, fyv);

    // printing values
    printf("%d %g %g\n", i, fx[i], fy[i]);

    // there is a slight difference in values of the original computation because of _mm256_rsqrt_ps. 
    // (I have used the same as mentioned in the hints) 
    // To obtain same result one could use _mm256_sqrt_ps and _mm256_div_ps of overall r^3,
    //  as the errors due to reciprocating and mutliplication is reduced
  }
}
