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
        let log2n = logn % 2;
        let log4n = logn / 2;

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

        Self {
            n,
            logn,
            log2n,
            log4n,
            ws,
            qw,
        }
    }
}

impl super::FFT for Radix4 {
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
        for _ in 0..self.log4n {
            m /= 4;
            for i in 0..m {
                let w1 = &self.ws[iw + i];
                let w2 = &self.ws[iw + 2 * m + i];
                let w3 = &(*w1 * w2);
                for j in 0..l {
                    let k0 = j * 4 * m + i;
                    let k1 = j * 4 * m + i + m;
                    let k2 = j * 4 * m + i + 2 * m;
                    let k3 = j * 4 * m + i + 3 * m;
                    let x0 = &x[k0];
                    let x2 = &x[k2];
                    let x1 = &x[k1];
                    let x3 = &x[k3];
                    let y0 = *x0 + x2;
                    let y1 = *x1 + x3;
                    let y2 = *x0 - x2;
                    let y3 = (*x1 - x3).i();
                    x[k0] = y0 + &y1;
                    x[k1] = (y0 - &y1) * w2;
                    x[k2] = (y2 + &y3) * w1;
                    x[k3] = (y2 - &y3) * w3;
                }
            }
            l *= 4;
            iw += 3 * m;
        }

        for _ in 0..self.log2n {
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
