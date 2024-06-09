use std::f64::consts::PI;

use super::Complex;

pub struct DITFactory {}

impl super::FFTFactory for DITFactory {
    fn create(&self, logn: usize) -> Box<dyn super::FFT> {
        Box::new(DIT::new(logn))
    }

    fn name(&self) -> &str {
        "DIT"
    }
}

struct DIT {
    n: usize,
    logn: usize,
    log2n: usize,
    log4n: usize,
}

impl DIT {
    pub fn new(logn: usize) -> Self {
        Self {
            n: 1 << logn,
            logn,
            log2n: logn % 2,
            log4n: logn / 2,
        }
    }
}

impl super::FFT for DIT {
    fn rft(&self, x: &mut Vec<Complex>) {
        let t = 2.0 * PI / self.n as f64 / 4.0;
        for (i, xi) in x.iter_mut().enumerate() {
            let theta = t * i as f64;
            let w = Complex::from((theta.cos(), theta.sin()));
            *xi *= &w;
        }
        self.dft(x);
    }

    fn irft(&self, x: &mut Vec<Complex>) {
        self.idft(x);

        let t = 2.0 * PI / self.n as f64 / 4.0;
        for (i, xi) in x.iter_mut().enumerate() {
            let theta = t * i as f64;
            let w = Complex::from((theta.cos(), theta.sin())).conj();
            *xi *= &w;
        }
    }

    fn dft(&self, x: &mut Vec<Complex>) {
        let mut m = self.n;
        let mut l = 1;
        let mut theta = 2.0 * PI / self.n as f64;
        for _ in 0..self.logn {
            m /= 2;
            for k in 0..l {
                for j in 0..m {
                    let t = theta * j as f64;
                    let w1 = Complex::from((t.cos(), t.sin()));
                    let j0 = 2 * m * k + j;
                    let j1 = 2 * m * k + j + m;
                    let x0 = x[j0].clone();
                    let x1 = x[j1].clone();
                    x[j0] = x0 + &x1;
                    x[j1] = (x0 - &x1) * &w1;
                }
            }
            theta *= 2.0;
            l *= 2;
        }
    }

    fn idft(&self, x: &mut Vec<Complex>) {
        let mut m = 1;
        let mut l = self.n;
        let mut theta = 2.0 * PI;
        for _ in 0..self.logn {
            l /= 2;
            theta *= 0.5;
            for k in 0..l {
                for j in 0..m {
                    let t = theta * j as f64;
                    let w1 = Complex::from((t.cos(), t.sin())).conj();
                    let j0 = 2 * m * k + j;
                    let j1 = 2 * m * k + j + m;
                    let x0 = x[j0].clone();
                    let x1 = x[j1].clone();
                    x[j0] = x0 + &x1;
                    x[j1] = (x0 - &x1) * &w1;
                }
            }
            m *= 2;
        }
    }

    fn bytes(&self) -> usize {
        0
    }

    fn flops(&self) -> usize {
        self.n / 2 * self.logn * 10
    }
}