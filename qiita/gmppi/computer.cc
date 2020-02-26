#include "computer.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

std::unique_ptr<Computer> Computer::Create(const Configuration& config) {
  switch (config.formula) {
    case Formula::kChudnovsky:
      return std::make_unique<Chudnovsky>(config);
    case Formula::kRamanujan:
      return std::make_unique<Ramanujan>(config);
  }
  return std::make_unique<Chudnovsky>(config);
}

void Computer::drm(const int64_t n0,
                   const int64_t n1,
                   Parameter& p0,
                   bool need_z) {
  if (n0 + 1 == n1) {
    setXYZ(n0, p0);
    return;
  }

  int64_t m = (n0 + n1) / 2;
  Parameter p1;
  drm(n0, m, p0, true);
  drm(m, n1, p1, need_z);

  // y0 = x1 * y0 + y1 * z0;
  p0.y *= p1.x;
  p1.y *= p0.z;
  p0.y += p1.y;

  p0.x *= p1.x;
  if (need_z)
    p0.z *= p1.z;
}

void Computer::parallel_drm(const int64_t n0,
                            const int64_t n1,
                            Parameter& p0,
                            bool need_z,
                            const int number_of_threads) {
  if (number_of_threads == 1) {
    drm(n0, n1, p0, need_z);
    return;
  }

  int num_new_threads = number_of_threads / 2;
  int num_rest_threads = number_of_threads - num_new_threads;

  int64_t m = (n0 + n1) / 2;
  Parameter p1;
  {
    auto clausure = [&] { parallel_drm(m, n1, p1, need_z, num_new_threads); };
    std::thread th(clausure);
    parallel_drm(n0, m, p0, true, num_rest_threads);
    th.join();
  }

  // y0 = x1 * y0 + y1 * z0;
  p0.y *= p1.x;
  p1.y *= p0.z;
  p0.y += p1.y;

  p0.x *= p1.x;
  if (need_z)
    p0.z *= p1.z;
}

void Computer::compute() {
  const int64_t precision = config_.digits * std::log2(10) + 10;
  pi_.set_prec(precision);

  const int64_t n = terms(config_.digits);
  Parameter param;
  parallel_drm(0, n, param, false, config_.number_of_threads);
  postProcess(param);
}

void Computer::output() {
  if (!config_.outputs)
    return;

  char filename[50];
  std::snprintf(filename, 45, "pi_peria_%s.out", name());
  std::ofstream ofs(filename);
  if (!ofs.is_open()) {
    std::cerr << "Could not open the file\n";
    return;
  }
  ofs << std::setprecision(config_.digits + 2) << pi_;
}

void Computer::check() {
  if (!config_.compare_file.has_value())
    return;

  const int64_t precision = config_.digits * std::log2(10) + 10;
  std::ifstream ifs(config_.compare_file.value());
  if (!ifs.is_open())
    return;
  mpf_class answer(0, precision);
  ifs >> answer;
  answer -= pi_;

  static constexpr int kBase = 10;
  static constexpr int kDigits = 10;
  mpf_out_str(stderr, kBase, kDigits, answer.get_mpf_t());
}
