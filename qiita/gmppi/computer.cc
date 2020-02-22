#include "computer.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

void Computer::drm(const int64_t n0,
                   const int64_t n1,
                   mpz_class& x0,
                   mpz_class& y0,
                   mpz_class& z0,
                   bool need_z) {
  if (n0 + 1 == n1) {
    setXYZ(n0, x0, y0, z0);
    return;
  }

  int64_t m = (n0 + n1) / 2;
  mpz_class x1, y1, z1;
  drm(n0, m, x0, y0, z0, true);
  drm(m, n1, x1, y1, z1, need_z);

  // y0 = x1 * y0 + y1 * z0;
  y0 *= x1;
  y1 *= z0;
  y0 += y1;

  x0 *= x1;
  if (need_z)
    z0 *= z1;
}

void Computer::parallel_drm(const int64_t n0,
                            const int64_t n1,
                            mpz_class& x0,
                            mpz_class& y0,
                            mpz_class& z0,
                            bool need_z,
                            const int number_of_threads) {
  if (number_of_threads == 1) {
    drm(n0, n1, x0, y0, z0, need_z);
    return;
  }

  int num_new_threads = number_of_threads / 2;
  int num_rest_threads = number_of_threads - num_new_threads;

  int64_t m = (n0 + n1) / 2;
  mpz_class x1, y1, z1;
  {
    auto clausure = [&] {
      parallel_drm(m, n1, x1, y1, z1, need_z, num_new_threads);
    };
    std::thread th(clausure);
    parallel_drm(n0, m, x0, y0, z0, true, num_rest_threads);
    th.join();
  }

  // y0 = x1 * y0 + y1 * z0;
  y0 *= x1;
  y1 *= z0;
  y0 += y1;

  x0 *= x1;
  if (need_z)
    z0 *= z1;
}

void Computer::compute(int number_of_threads) {
  const int64_t precision = digits_ * std::log2(10) + 10;
  pi_.set_prec(precision);

  const int64_t n = terms(digits_);
  mpz_class x, y, z;
  parallel_drm(0, n, x, y, z, false, number_of_threads);
  postProcess(x, y);
}

void Computer::output() {
  char filename[50];
  std::snprintf(filename, 45, "pi_peria_%s.out", name());
  std::ofstream ofs(filename);
  if (!ofs.is_open()) {
    std::cerr << "Could not open the file\n";
    return;
  }
  ofs << std::setprecision(digits_ + 2) << pi_;
}

void Computer::check(const char* answer_filename) {
  if (!answer_filename)
    return;

  const int64_t precision = digits_ * std::log2(10) + 10;
  std::ifstream ifs(answer_filename);
  if (!ifs.is_open())
    return;
  mpf_class answer(0, precision);
  ifs >> answer;
  answer -= pi_;

  static constexpr int kBase = 10;
  static constexpr int kDigits = 10;
  mpf_out_str(stderr, kBase, kDigits, answer.get_mpf_t());
}
