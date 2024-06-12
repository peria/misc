mod complex;
mod radix2;
mod radix4;

pub use complex::Complex;
pub use radix2::Radix2Factory;
pub use radix4::Radix4Factory;

pub trait FFTFactory {
    fn create(&self, n: usize) -> Box<dyn FFT>;
    fn name(&self) -> &str;
}

pub trait FFT {
    fn rft(&self, x: &mut Vec<Complex>);
    fn irft(&self, x: &mut Vec<Complex>);
    fn dft(&self, x: &mut Vec<Complex>);
    fn idft(&self, x: &mut Vec<Complex>);
    fn bytes(&self) -> usize;
    fn flops(&self) -> usize;
}

#[cfg(test)]
mod test {
    use super::Complex;
    use super::FFTFactory;
    use super::Radix2Factory;
    use super::Radix4Factory;

    #[test]
    fn back_test() {
        let factories: Vec<Box<dyn FFTFactory>> =
            vec![Box::new(Radix2Factory {}), Box::new(Radix4Factory {})];
        for factory in factories.iter() {
            for logn in 2..10 {
                let n = 1 << logn;
                let mut x = vec![Complex::new(); n];
                for (i, xi) in x.iter_mut().enumerate() {
                    xi.real = i as f64 + 1.0;
                }
                let fft = factory.create(logn);
                fft.rft(&mut x);
                fft.irft(&mut x);
                for (i, xi) in x.iter().enumerate() {
                    eprintln!("[{:4}/{:4}] {:.4} - {:.4}", i, n, xi.real, xi.imag);
                    debug_assert!((xi.real - i as f64 - 1.0).abs() < 1e-3);
                    debug_assert!(xi.imag.abs() < 1e-3);
                }
            }
            eprintln!("{} OK", factory.name());
        }
    }

    #[test]
    fn convolution_test() {
        let factories: Vec<Box<dyn FFTFactory>> =
            vec![Box::new(Radix2Factory {}), Box::new(Radix4Factory {})];
        for factory in factories.iter() {
            eprintln!("Testing {}", factory.name());
            for logn in 2..10 {
                let n = 1 << logn;
                let mut x = vec![Complex::new(); n];
                for xi in x.iter_mut() {
                    xi.real = 1.0;
                }
                let fft = factory.create(logn);
                fft.rft(&mut x);
                for xi in x.iter_mut() {
                    let r = xi.real;
                    let i = xi.imag;
                    xi.real = r * r - i * i;
                    xi.imag = 2.0 * r * i;
                }
                fft.irft(&mut x);
                for (i, xi) in x.iter().enumerate() {
                    eprintln!("[{:4}/{:4}] {:.4} + {:.4}i", i, n, xi.real, xi.imag);
                }
                for (i, xi) in x.iter().enumerate() {
                    debug_assert!((xi.real - (1 + i) as f64).abs() < 1e-3);
                    debug_assert!((xi.imag - (n - 1 - i) as f64).abs() < 1e-3);
                }
            }
            eprintln!("{} OK", factory.name());
        }
    }
}
