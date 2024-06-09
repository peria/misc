mod complex;
mod radix2;

pub use complex::Complex;
pub use radix2::Radix2Factory;

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
