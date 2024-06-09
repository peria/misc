use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

#[derive(Copy, Clone)]
pub struct Complex {
    real: f64,
    imag: f64,
}

impl Complex {
    pub fn new() -> Self {
        Complex {
            real: 0.0,
            imag: 0.0,
        }
    }

    pub fn conj(&self) -> Self {
        Complex {
            real: self.real,
            imag: -self.imag,
        }
    }

    pub fn i(&self) -> Self {
        Complex {
            real: -self.imag,
            imag: self.real,
        }
    }
}

impl From<(f64, f64)> for Complex {
    fn from(value: (f64, f64)) -> Self {
        Complex {
            real: value.0,
            imag: value.1,
        }
    }
}

impl AddAssign<&Self> for Complex {
    fn add_assign(&mut self, rhs: &Self) {
        self.real += rhs.real;
        self.imag += rhs.imag;
    }
}

impl SubAssign<&Self> for Complex {
    fn sub_assign(&mut self, rhs: &Self) {
        self.real -= rhs.real;
        self.imag -= rhs.imag;
    }
}

impl MulAssign<&Self> for Complex {
    fn mul_assign(&mut self, rhs: &Self) {
        let real = self.real * rhs.real - self.imag * rhs.imag;
        let imag = self.imag * rhs.real + self.real * rhs.imag;
        self.real = real;
        self.imag = imag;
    }
}

impl MulAssign<f64> for Complex {
    fn mul_assign(&mut self, rhs: f64) {
        self.real *= rhs;
        self.imag *= rhs;
    }
}

impl Add<&Self> for Complex {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self::Output {
        Complex {
            real: self.real + rhs.real,
            imag: self.imag + rhs.imag,
        }
    }
}

impl Sub<&Self> for Complex {
    type Output = Self;
    fn sub(self, rhs: &Self) -> Self::Output {
        Complex {
            real: self.real - rhs.real,
            imag: self.imag - rhs.imag,
        }
    }
}

impl Mul<&Self> for Complex {
    type Output = Self;
    fn mul(self, rhs: &Self) -> Self::Output {
        Complex {
            real: self.real * rhs.real - self.imag * rhs.imag,
            imag: self.imag * rhs.real + self.real * rhs.imag,
        }
    }
}
