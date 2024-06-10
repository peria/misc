use std::f64::consts::PI;

use super::Complex;

pub struct Radix2Factory {}

impl super::FFTFactory for Radix2Factory {
    fn create(&self, logn: usize) -> Box<dyn super::FFT> {
        Box::new(Radix2::new(logn))
    }

    fn name(&self) -> &str {
        "R2"
    }
}

struct Radix2 {
    n: usize,
    logn: usize,
    ws: Vec<Complex>,
    qw: Complex,
}

impl Radix2 {
    pub fn new(logn: usize) -> Self {
        let n = 1 << logn;

        let mut ws = Vec::new();
        let mut m = n;
        let mut theta = 2.0 * PI / n as f64;
        for _ in 0..logn {
            m /= 2;
            for i in 0..m {
                let t = theta * i as f64;
                ws.push(Complex::from((t.cos(), t.sin())));
            }
            theta *= 2.0;
        }

        let qt = 2.0 * PI / n as f64 / 4.0;
        let qw = Complex::from((qt.cos(), qt.sin()));

        Self { n, logn, ws, qw }
    }
}

impl super::FFT for Radix2 {
    fn rft(&self, x: &mut Vec<Complex>) {
        let qw1 = &self.qw;
        let qw2 = &(*qw1 * qw1);
        let qw3 = &(*qw2 * qw1);
        for i in 0..(self.n / 4) {
            let w0 = &self.ws[i];
            let w1 = &(*w0 * qw1);
            let w2 = &(*w0 * qw2);
            let w3 = &(*w0 * qw3);
            x[4 * i + 0] *= w0;
            x[4 * i + 1] *= w1;
            x[4 * i + 2] *= w2;
            x[4 * i + 3] *= w3;
        }
        self.dft(x);
    }

    fn irft(&self, x: &mut Vec<Complex>) {
        self.idft(x);

        let qw1 = &self.qw.conj();
        let qw2 = &(*qw1 * qw1);
        let qw3 = &(*qw2 * qw1);
        for i in 0..(self.n / 4) {
            let w0 = &self.ws[i].conj();
            let w1 = &(*w0 * qw1);
            let w2 = &(*w0 * qw2);
            let w3 = &(*w0 * qw3);
            x[4 * i + 0] *= w0;
            x[4 * i + 1] *= w1;
            x[4 * i + 2] *= w2;
            x[4 * i + 3] *= w3;
        }
    }

    fn dft(&self, x: &mut Vec<Complex>) {
        let mut m = self.n;
        let mut l = 1;
        let mut iw = 0;
        for _ in 0..self.logn {
            m /= 2;
            for i in 0..m {
                let w1 = &self.ws[iw + i];
                for j in 0..l {
                    let k0 = j * 2 * m + i;
                    let k1 = j * 2 * m + i + m;
                    let x0 = x[k0].clone();
                    let x1 = x[k1].clone();
                    x[k0] = x0 + &x1;
                    x[k1] = (x0 - &x1) * w1;
                }
            }
            l *= 2;
            iw += m;
        }
    }

    fn idft(&self, x: &mut Vec<Complex>) {
        let mut m = 1;
        let mut l = self.n;
        let mut iw = self.n - 1;
        for _ in 0..self.logn {
            l /= 2;
            iw -= m;
            for i in 0..m {
                let w1 = &self.ws[iw + i].conj();
                for j in 0..l {
                    let k0 = j * 2 * m + i;
                    let k1 = j * 2 * m + i + m;
                    let x0 = x[k0].clone();
                    let x1 = x[k1] * w1;
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
