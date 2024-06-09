use std::f64::consts::PI;

use super::Complex;

pub struct Radix4Factory {}

impl super::FFTFactory for Radix4Factory {
    fn create(&self, logn: usize) -> Box<dyn super::FFT> {
        Box::new(Radix4::new(logn))
    }

    fn name(&self) -> &str {
        "R4"
    }
}

struct Radix4 {
    n: usize,
    logn: usize,
    log2n: usize,
    log4n: usize,
    ws: Vec<Complex>,
    qw: Complex,
}

impl Radix4 {
    pub fn new(logn: usize) -> Self {
        let n = 1 << logn;

        let mut ws = Vec::new();
        let mut m = n;
        let mut theta = 2.0 * PI / n as f64;
        for _ in 0..logn {
            m /= 2;
            for j in 0..m {
                let t = theta * j as f64;
                let w1 = Complex::from((t.cos(), t.sin()));
                ws.push(w1);
            }
            theta *= 2.0;
        }
        let qt = theta / 4.0;
        let qw = Complex::from((qt.cos(), qt.sin()));

        Self {
            n,
            logn,
            log2n: logn % 2,
            log4n: logn / 2,
            ws,
            qw,
        }
    }
}

impl super::FFT for Radix4 {
    fn rft(&self, x: &mut Vec<Complex>) {
        for i in 0..(self.n / 4) {
            let w0 = &self.ws[i];
            let w1 = &(self.qw * w0);
            let w2 = &(self.qw * w1);
            let w3 = &(self.qw * w2);
            x[4 * i + 0] *= w0;
            x[4 * i + 1] *= w1;
            x[4 * i + 2] *= w2;
            x[4 * i + 3] *= w3;
        }
        self.dft(x);
    }

    fn irft(&self, x: &mut Vec<Complex>) {
        self.idft(x);

        let qw = self.qw.conj();
        for i in 0..(self.n / 4) {
            let w0 = &self.ws[i].conj();
            let w1 = &(qw * w0);
            let w2 = &(qw * w1);
            let w3 = &(qw * w2);
            x[4 * i + 0] *= w0;
            x[4 * i + 1] *= w1;
            x[4 * i + 2] *= w2;
            x[4 * i + 3] *= w3;
        }
    }

    fn dft(&self, x: &mut Vec<Complex>) {
        let mut m = self.n;
        let mut l = 1;
        let mut theta = 2.0 * PI / self.n as f64;
        for _ in 0..self.logn {
            m /= 2;
            for i in 0..l {
                let t = theta * i as f64;
                let w1 = Complex::from((t.cos(), t.sin()));
                for j in 0..m {
                    let k0 = j + i * 2 * m;
                    let k1 = i * 2 * m + j + m;
                    let x0 = x[k0].clone();
                    let x1 = x[k1].clone();
                    x[k0] = x0 + &x1;
                    x[k1] = (x0 - &x1) * &w1;
                }
            }
            l *= 2;
            theta *= 2.0;
        }
    }

    fn idft(&self, x: &mut Vec<Complex>) {
        let mut m = 1;
        let mut l = self.n;
        let mut theta = 2.0 * PI;
        for _ in 0..self.logn {
            l /= 2;
            theta /= 2.0;
            for i in 0..l {
                let t = theta * i as f64;
                let w1 = Complex::from((t.cos(), t.sin())).conj();
                for j in 0..m {
                    let k0 = j + i * 2 * m;
                    let k1 = i * 2 * m + j + m;
                    let x0 = x[k0].clone();
                    let x1 = x[k1] * &w1;
                    x[k0] = x0 + &x1;
                    x[k1] = x0 - &x1;
                }
            }
            m *= 2;
        }

        let inv = 1.0 / self.n as f64;
        x.iter_mut().for_each(|xi| *xi *= inv);
    }

    fn bytes(&self) -> usize {
        self.ws.len() * 16 + 16
    }

    fn flops(&self) -> usize {
        self.n / 2 * self.logn * 10
    }
}
