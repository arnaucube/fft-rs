# fft-rs [![Test](https://github.com/arnaucube/fft-rs/workflows/Test/badge.svg)](https://github.com/arnaucube/fft-rs/actions?query=workflow%3ATest)

Fast Fourier Transform implementation in Rust.

https://en.wikipedia.org/wiki/Fast_Fourier_transform & [DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)

## Usage
```rust
let values: Vec<f64> = vec![0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];

// compute the FFT (Fast Fourier Transform)
let fft_res = fft(&values);

// compute the IFFT (Inverse Fast Fourier Transform)
let ifft_res = ifft(&fft_res);


// Also, available directly (and slower) DFT & IDFT:

// compute the DFT (Discrete Fourier Transform)
let dft_res = dft(&values);

// compute the IDFT (Inverse Discrete Fourier Transform)
let idft_res = idft(&dft_res);
```
