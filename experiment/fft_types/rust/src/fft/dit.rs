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
    ws: Vec<Complex>,
    qw: Complex,
}

impl DIT {
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

impl super::FFT for DIT {
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
        let mut wi = 0;
        for _ in 0..self.logn {
            m /= 2;
            for k in 0..l {
                for j in 0..m {
                    let w1 = &self.ws[wi + j];
                    let j0 = 2 * m * k + j;
                    let j1 = 2 * m * k + j + m;
                    let x0 = x[j0].clone();
                    x[j0] = x0 + &x[j1];
                    x[j1] = (x0 - &x[j1]) * w1;
                }
            }
            wi += m;
            l *= 2;
        }
    }

    fn idft(&self, x: &mut Vec<Complex>) {
        let mut m = 1;
        let mut l = self.n;
        let mut wi = self.ws.len();
        for _ in 0..self.logn {
            l /= 2;
            wi -= m;
            for k in 0..l {
                for j in 0..m {
                    let w1 = self.ws[wi + j].conj();
                    let j0 = 2 * m * k + j;
                    let j1 = 2 * m * k + j + m;
                    let x0 = x[j0].clone();
                    x[j0] = x0 + &x[j1];
                    x[j1] = (x0 - &x[j1]) * &w1;
                }
            }
            m *= 2;
        }
    }

    fn bytes(&self) -> usize {
        self.ws.len() * 16 + 16
    }

    fn flops(&self) -> usize {
        self.n / 2 * self.logn * 10
    }
}
