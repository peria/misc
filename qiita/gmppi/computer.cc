#include "computer.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

Computer::Factorizer::Factorizer(const int64_t n) : elements_(n / 2 + 1) {
  elements_[1 / 2].prime_factor = PrimeFactor{1, 1};
  for (int64_t p = 3; p <= n; p += 2) {
    Element& prime = elements_[p / 2];
    if (prime.prime_factor.prime)
      continue;
    prime.prime_factor = PrimeFactor{p, 1};
    for (int64_t c = p * p, next_index = p / 2; c <= n;
         c += p * 2, ++next_index) {
      Element& composite = elements_[c / 2];
      if (composite.prime_factor.prime)
        continue;
      composite.prime_factor.prime = p;
      Element& next = elements_[next_index];
      if (next.prime_factor.prime == p) {
        composite.prime_factor.exp = next.prime_factor.exp + 1;
        composite.next = next.next;
      } else {
        composite.prime_factor.exp = 1;
        composite.next = next_index;
      }
    }
  }
}

Computer::Factorized Computer::Factorizer::factorize(int64_t n, int64_t exp) {
  Factorized factors;
  while ((n & 1) == 0)
    n >>= 1;
  for (int64_t i = n / 2; i; i = elements_[i].next) {
    auto& pf = elements_[i].prime_factor;
    factors.add(pf);
  }
  return factors;
}

Computer::Factorized& Computer::Factorized::operator*=(
    const Factorized& other) {
  std::vector<PrimeFactor> buffer;
  std::swap(buffer, factors);
  auto i = buffer.begin();
  auto j = other.factors.begin();
  while (i != buffer.end() && j != other.factors.end()) {
    if (i->prime < j->prime) {
      add(*i);
      ++i;
    } else if (i->prime > j->prime) {
      add(*j);
      ++j;
    } else {
      add(PrimeFactor{i->prime, i->exp + j->exp});
      ++i;
      ++j;
    }
  }
  for (; i != buffer.end(); ++i)
    add(*i);
  for (; j != other.factors.end(); ++j)
    add(*j);
  return *this;
}

mpz_class Computer::Factorized::toMpz() const {
  return toMpz(0, factors.size());
}

mpz_class Computer::Factorized::toMpz(int64_t a, int64_t b) const {
  if (b - a <= 3) {
    mpz_class x = 1;
    for (int64_t i = a; i < b; ++i) {
      for (int64_t j = 0; j < factors[i].exp; ++j) {
        x *= factors[i].prime;
      }
    }
    return x;
  }
  int64_t m = (a + b) / 2;
  return toMpz(a, m) * toMpz(m, b);
}

void Computer::Factorized::shrink() {
  auto i = factors.begin();
  for (; i != factors.end() && i->exp; ++i) {
  }

  auto j = i;
  for (; i != factors.end(); ++i) {
    if (i->exp)
      *j++ = *i;
  }
  factors.erase(j, factors.end());
}

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

  reduceGcd(p0, p1);

  // y0 = x1 * y0 + y1 * z0;
  p0.y *= p1.x;
  p1.y *= p0.z;
  p0.y += p1.y;

  p0.x *= p1.x;
  if (need_z)
    p0.z *= p1.z;
}

void Computer::parallelDrm(const int64_t n0,
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
    auto clausure = [&] { parallelDrm(m, n1, p1, need_z, num_new_threads); };
    std::thread th(clausure);
    parallelDrm(n0, m, p0, true, num_rest_threads);
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

void Computer::reduceGcd(Parameter& p0, Parameter& p1) {
  Factorized f_gcd = splitGcd(p0.factorized_z, p1.factorized_x);
  mpz_class g = f_gcd.toMpz();
  mpz_divexact(p0.z.get_mpz_t(), p0.z.get_mpz_t(), g.get_mpz_t());
  mpz_divexact(p1.x.get_mpz_t(), p1.x.get_mpz_t(), g.get_mpz_t());
}

Computer::Factorized Computer::splitGcd(Factorized& a, Factorized& b) {
  Factorized g;
  for (auto ia = a.factors.begin(), ib = b.factors.begin();
       ia != a.factors.end() && ib != b.factors.end();) {
    if (ia->prime < ib->prime) {
      ++ia;
    } else if (ia->prime > ib->prime) {
      ++ib;
    } else {
      int64_t e = std::min(ia->exp, ib->exp);
      ia->exp -= e;
      ib->exp -= e;
      g.factors.push_back(PrimeFactor{ia->prime, e});
      ++ia;
      ++ib;
    }
  }
  a.shrink();
  b.shrink();
  return g;
}

void Computer::compute() {
  const int64_t precision = config_.digits * std::log2(10) + 10;
  pi_.set_prec(precision);

  const int64_t n = terms();
  Parameter param;
  parallelDrm(0, n, param, false, config_.number_of_threads);
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
