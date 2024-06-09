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

mod test {
    use super::Complex;
    use super::FFTFactory;
    use super::Radix2Factory;
    use super::Radix4Factory;

    #[test]
    fn reverse() {
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
}
