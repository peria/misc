#include "fft.h"

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

namespace {

double* g_work;
double* g_table;

}  // namespace

void RFT::SetGlobals(void* work, void* table, void*) {
  g_work = static_cast<double*>(work);
  g_table = static_cast<double*>(table);
}

bool RFT::Validate(double* data) {
  for (int log2n = 1; log2n < 10; ++log2n) {
    const int n = 1 << log2n;
    InitTable(log2n);
    for (int i = 0; i < n; ++i)
      data[i] = i;
    Forward(log2n, n, data);
    Backward(log2n, n, data);
    for (int i = 0; i < n; ++i) {
      if (abs(data[i] - i) > 1e-4) {
        cout << "2^" << log2n << " [" << i << "] " << data[i] << "\n";
        return false;
      }
    }
  }
  return true;
}

void FFT::InitTable(const int log2n) {
  double* ptr = g_table;

  int width = 1;
  for (int k = 0; k < log2n; ++k) {
    for (int i = 0; i < width; ++i) {
      double wt = M_PI * i / width;
      ptr[2*i] = cos(wt);
      ptr[2*i+1] = -sin(wt);
    }
    ptr += 2 * width;
    width *= 2;
  }
}

void RFT::InitTable(const int log2n) {
  FFT::InitTable(log2n - 1);
}

void FFT::Radix2(const int width, const int height,
                 double* ptr, double* x, double* y) {
  for (int i = 0; i < width; ++i) {
    double wr = ptr[2*i], wi = ptr[2*i+1];
    for (int j = 0; j < height; ++j) {
      int ix0 = j * width + i, ix1 = ix0 + height * width;
      int iy0 = j * 2 * width + i, iy1 = iy0 + width;
      double tr = x[2*ix1] * wr - x[2*ix1+1] * wi;
      double ti = x[2*ix1] * wi + x[2*ix1+1] * wr;
      y[2*iy1  ] = x[2*ix0  ] - tr;
      y[2*iy1+1] = x[2*ix0+1] - ti;
      y[2*iy0  ] = x[2*ix0  ] + tr;
      y[2*iy0+1] = x[2*ix0+1] + ti;
    }
  }
}

void FFT::Forward(const int log2n, const int n, double* data) {
  double* x = data;
  double* y = g_work;
  double* table = g_table;

  int width = 1, height = n;
  for (int i = 0; i < log2n; ++i) {
    height /= 2;
    if (i % 2) {
      Radix2(width, height, table, y, x);
    } else {
      Radix2(width, height, table, x, (i == log2n - 1) ? x : y);
    }
    table += 2 * width;
    width *= 2;
  }
}

void FFT::Backward(const int log2n, const int n, double* data) {
  for (int i = 0; i < n; ++i) {
    data[2*i+1] = -data[2*i+1];
  }
  Forward(log2n, n, data);
  double inv = 1.0 / n;
  for (int i = 0; i < n; ++i) {
    data[2*i] *= inv;
    data[2*i+1] *= -inv;
  }
}

void RFT::Forward(const int log2n, const int real_n, double* x) {
  const int n = real_n / 2;

  FFT::Forward(log2n - 1, n, x);
  // Convert from complex to real
  double t = 2 * M_PI / real_n;
  for (int k = 1, j = n - 1; k < n / 2; ++k, --j) {
    double x0r = x[2*k], x0i = x[2*k+1];
    double x1r = x[2*j], x1i = x[2*j+1];
    double ar = x0r - x1r, ai = x0i + x1i;
    double wr = 1 + sin(k * t);
    double wi = cos(k * t);
    double br = (ar * wr - ai * wi) * 0.5;
    double bi = (ar * wi + ai * wr) * 0.5;
    x[2*k] = x0r - br;
    x[2*k+1] = x0i - bi;
    x[2*j] = x1r + br;
    x[2*j+1] = -(-x1i + bi);
  }
  double xr = x[0], xi = x[1];
  x[0] = xr + xi;
  x[1] = xr - xi;
}

void RFT::Backward(const int log2n, const int real_n, double* x) {
  const int n = real_n / 2;

  // Convert from real to complex
  double t = M_PI / n;
  for (int k = 1, j = n - 1; k < n / 2; ++k, --j) {
    double x0r = x[2*k], x0i = x[2*k+1];
    double x1r = x[2*j], x1i = x[2*j+1];
    double ar = x0r - x1r, ai = x0i + x1i;
    double wr = 1 + sin(k * t);
    double wi = -cos(k * t);
    double br = (ar * wr - ai * wi) * 0.5;
    double bi = (ar * wi + ai * wr) * 0.5;
    x[2*k] = x0r - br;
    x[2*k+1] = x0i - bi;
    x[2*j] = x1r + br;
    x[2*j+1] = -(-x1i + bi);
  }
  double xr = x[0], xi = x[1];
  x[0] = (xr + xi) * 0.5;
  x[1] = (xr - xi) * 0.5;

  FFT::Backward(log2n - 1, n, x);
}
