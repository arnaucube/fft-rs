use criterion::{criterion_group, criterion_main, Criterion};

use num::complex::{Complex, Complex64};
extern crate rand;

use rand::Rng;

use fft_rs;

fn criterion_benchmark(c: &mut Criterion) {
    let values: Vec<f64> = rand::thread_rng()
        .sample_iter(rand::distributions::Standard)
        .take(1024)
        .collect();

    c.bench_function("dft", |b| b.iter(|| fft_rs::dft(&values)));

    c.bench_function("fft", |b| b.iter(|| fft_rs::fft(&values)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
