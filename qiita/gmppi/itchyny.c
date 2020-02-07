/* Calculate pi based on Chudnovsky algorithm (& Binary Splitting method), using GMP */
/* [1] Computation of 2700 billion decimal digits of Pi using a Desktop Computer,
       Fabrice Bellard, Feb 11 2010 (4th revision),
       http://bellard.org/pi/pi2700e9/pipcrecord.pdf */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>
#include <time.h>

#define DEFAULT_DIGITS 10000

mpz_t A, B, C, C3over24;

void computePQT(mpz_t P, mpz_t Q, mpz_t T, int n1, int n2) {
  int m;
  if (n1 + 1 == n2) {
    mpz_t a;
    mpz_init(a);
    /* P(n2-1, n2) = - (2 * n2 - 1) * (6 * n2 - 5) * (6 * n2 - 1) */
    mpz_set_ui(P, (2 * n2 - 1));
    mpz_mul_ui(P, P, (6 * n2 - 5));
    mpz_mul_ui(P, P, (6 * n2 - 1));
    mpz_neg(P, P);
    /* Q(n2-1, n2) = C^3 * n2^3 / 24 */
    mpz_mul_ui(Q, C3over24, n2);
    mpz_mul_ui(Q, Q, n2);
    mpz_mul_ui(Q, Q, n2);
    /* a_n2 = (A + B * n2) */
    mpz_mul_ui(a, B, n2);
    mpz_add(a, a, A);
    /* T(n2-1, n2) = a_n2 * P(n2-1, n2) */
    mpz_mul(T, a, P);
    mpz_clear(a);
  } else {
    mpz_t P1, Q1, T1, P2, Q2, T2;
    mpz_init(P1);
    mpz_init(Q1);
    mpz_init(T1);
    mpz_init(P2);
    mpz_init(Q2);
    mpz_init(T2);
    m = (n1 + n2) / 2;
    computePQT(P1, Q1, T1, n1, m);
    computePQT(P2, Q2, T2, m, n2);
    /* P = P1 * P2 */
    mpz_mul(P, P1, P2);
    /* Q = Q1 * Q2 */
    mpz_mul(Q, Q1, Q2);
    /* T = T1 * Q2 + P1 * T2 */
    mpz_mul(T, T1, Q2);
    mpz_mul(P1, P1, T2);
    mpz_add(T, T, P1);
    mpz_clear(P1);
    mpz_clear(Q1);
    mpz_clear(T1);
    mpz_clear(P2);
    mpz_clear(Q2);
    mpz_clear(T2);
  }
}

int main (int argc, char* argv[]) {
  clock_t start, end;
  FILE* fp;

  int digits = 0;
  if (argc > 1)
    digits = atoi(argv[1]);
  if (digits <= 0)
    digits = DEFAULT_DIGITS;

  int prec = digits * log2(10);
  int n = digits / 14;

  start = clock();

  mpf_t pi, temp;
  mpz_t P, Q, T;

  /* Initialize GMP numbers */
  mpf_set_default_prec(prec);
  mpf_init(pi);
  mpf_init(temp);
  mpz_init(P);
  mpz_init(Q);
  mpz_init(T);
  mpz_init(A);
  mpz_init(B);
  mpz_init(C);
  mpz_init(C3over24);

  /* Assignment */
  mpz_set_str(A, "13591409", 10);
  mpz_set_str(B, "545140134", 10);
  mpz_set_str(C, "640320", 10);
  mpz_mul(C3over24, C, C);
  mpz_mul(C3over24, C3over24, C);
  mpz_div_ui(C3over24, C3over24, 24);

  computePQT(P, Q, T, 0, n);

  /* pi = C ^ (1 / 2)     */
  mpf_set_z(temp, C);
  mpf_sqrt(pi, temp);
  /*      * C             */
  mpf_mul(pi, pi, temp);
  /*      * Q             */
  mpf_set_z(temp, Q);
  mpf_mul(pi, pi, temp);
  /*      / (T + A * Q)   */
  mpz_mul(Q, A, Q);
  mpz_add(Q, Q, T);
  mpf_set_z(temp, Q);
  mpf_div(pi, pi, temp);
  /*      / 12            */
  mpf_div_ui(pi, pi, 12);

  /* mpf_out_str(stdout, 10, digits, pi); */
  fp = fopen("pi_itchyny.out", "w");
  if (fp == NULL) {
    printf("Could not open the file.\n");
    exit(EXIT_FAILURE);
  }
  end = clock();
  printf("Compute: %.3f s\n",(double)(end - start) / CLOCKS_PER_SEC);
  mpf_out_str(fp, 10, digits, pi);

  mpf_clear(pi);
  mpf_clear(temp);
  mpz_clear(P);
  mpz_clear(Q);
  mpz_clear(T);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(C);
  mpz_clear(C3over24);

  end = clock();
  printf("All    : %.3f s\n",(double)(end - start) / CLOCKS_PER_SEC);

  return 0;
}
