mod fft;

use std::time;

use fft::Complex;
use fft::FFTFactory;

fn main() {
    let factories: Vec<Box<dyn FFTFactory>> = vec![
        Box::new(fft::Radix2Factory {}),
        Box::new(fft::Radix4Factory {}),
    ];

    print_names(&factories);
    measure_performace(&factories);
}

fn print_names(factories: &Vec<Box<dyn FFTFactory>>) {
    print!("     |");
    for factory in factories.iter() {
        print!(" {:5} |", factory.name());
    }
    println!("");
}

fn measure_performace(factories: &Vec<Box<dyn FFTFactory>>) {
    const MAX_LOGN: usize = 17;
    for logn in 2..=MAX_LOGN {
        print!("2^{:2} |", logn);
        let n = 1 << logn;
        let mut x = vec![Complex::new(); n];
        for factory in factories.iter() {
            let fft = factory.create(logn);
            let start = time::Instant::now();
            let due = start + time::Duration::from_secs_f64(2.0);
            let mut count = 0;
            while count < 10000 && time::Instant::now() < due {
                fft.rft(&mut x);
                fft.irft(&mut x);
                count += 1;
            }
            let duration = (time::Instant::now() - start).as_secs_f64();
            let mflops = 2.0 * fft.flops() as f64 * count as f64 / duration * 1e-6;
            print!("{:6} |", mflops.round() as i64);
        }
        println!("");
    }
}
