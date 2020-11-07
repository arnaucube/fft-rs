use num::complex::{Complex, Complex64};
use std::f64::consts::PI;

// fft computes the Fast Fourier Transform
pub fn fft(x: &Vec<f64>) -> Vec<Complex64> {
    let mut x_compl: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); x.len()];
    for i in 0..x.len() {
        x_compl[i] = Complex::new(x[i], 0_f64);
    }
    fft_compl(x_compl)
}

fn fft_compl(x: Vec<Complex64>) -> Vec<Complex64> {
    let N = x.len();
    if N % 2 > 0 {
        panic!("not a power of 2");
    } else if N <= 2 {
        return dft_compl(x);
    }

    let mut x_even: Vec<Complex64> = Vec::new();
    let mut x_odd: Vec<Complex64> = Vec::new();
    for i in 0..x.len() {
        if i % 2 == 0 {
            x_even.push(x[i]);
        } else {
            x_odd.push(x[i]);
        }
    }
    let mut x_even_cmplx = fft_compl(x_even);
    let mut x_odd_cmplx = fft_compl(x_odd);

    let mut w = Complex::new(0_f64, 2_f64 * PI / N as f64);
    let mut f_k: Vec<Complex64> = Vec::new();
    for k in 0..x.len() {
        let k_compl = Complex::new(k as f64, 0_f64);
        f_k.push((w * k_compl).exp());
    }

    let mut r: Vec<Complex64> = Vec::new();
    let mut aa = add_vv(
        x_even_cmplx.clone(),
        mul_vv_el(x_odd_cmplx.clone(), f_k.clone()[0..x.len() / 2].to_vec()),
    );
    let mut bb = add_vv(
        x_even_cmplx.clone(),
        mul_vv_el(x_odd_cmplx.clone(), f_k.clone()[x.len() / 2..].to_vec()),
    );
    r.append(&mut aa);
    r.append(&mut bb);

    r
}

// ifft computes the Inverse Fast Fourier Transform
pub fn ifft(x: &Vec<Complex64>) -> Vec<f64> {
    // use the IFFT method of computing conjugates, then FFT, then conjugate again, and then divide
    // by N
    let mut x_conj: Vec<Complex64> = Vec::new();
    for i in 0..x.len() {
        x_conj.push(x[i].conj());
    }

    let x_res = fft_compl(x_conj);

    let mut r: Vec<Complex64> = Vec::new();
    for i in 0..x_res.len() {
        r.push(x_res[i].conj());
    }

    let n = x.len() as f64;
    let mut rr: Vec<f64> = Vec::new();
    for i in 0..r.len() {
        r[i] = r[i] / Complex::new(n, 0_f64);
        rr.push(r[i].re);
    }
    rr
}

// dft computes the Discrete Fourier Transform
pub fn dft(x: &Vec<f64>) -> Vec<Complex64> {
    let mut x_compl: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); x.len()];
    for i in 0..x.len() {
        x_compl[i] = Complex::new(x[i], 0_f64);
    }
    dft_compl(x_compl)
}

fn dft_compl(x: Vec<Complex64>) -> Vec<Complex64> {
    let mut w = Complex::new(0_f64, -2_f64 * PI / x.len() as f64);

    // f_k (dft_matrix) = SUM{n=0, N-1} f_n * e^(-j2pi*k*n)/N
    // https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    let mut dft_matrix: Vec<Vec<Complex64>> = Vec::new();
    for i in 0..x.len() {
        let mut f_k: Vec<Complex64> = Vec::new();
        for j in 0..x.len() {
            let i_compl = Complex::new(0_f64, i as f64);
            let j_compl = Complex::new(0_f64, j as f64);
            let fe = (w * i_compl * j_compl).exp();
            f_k.push(fe);
        }
        dft_matrix.push(f_k.clone());
    }
    let r = mul_mv(dft_matrix, x);
    r
}

// idft computes the Inverse Discrete Fourier Transform
pub fn idft(x: &Vec<Complex64>) -> Vec<f64> {
    let mut w = Complex::new(0_f64, 2_f64 * PI / x.len() as f64);

    // f_k (dft_matrix) = (SUM{n=0, N-1} f_n * e^(j2pi*k*n)/N)/N
    let mut dft_matrix: Vec<Vec<Complex64>> = Vec::new();
    for i in 0..x.len() {
        let mut f_k: Vec<Complex64> = Vec::new();
        for j in 0..x.len() {
            let i_compl = Complex::new(0_f64, i as f64);
            let j_compl = Complex::new(0_f64, j as f64);
            let fe = (w * i_compl * j_compl).exp();
            f_k.push(fe);
        }
        dft_matrix.push(f_k.clone());
    }
    let mut r = mul_mv(dft_matrix, x.clone());
    let n = x.len() as f64;
    let mut rr: Vec<f64> = Vec::new();
    for i in 0..r.len() {
        r[i] = r[i] / Complex::new(n, 0_f64);
        rr.push(r[i].re);
    }
    rr
}

// mul_mv multiplies a Matrix by a Vector
fn mul_mv(a: Vec<Vec<Complex64>>, b: Vec<Complex64>) -> Vec<Complex64> {
    if a[0].len() != a.len() {
        panic!("err b[0].len():{:?} != b.len():{:?}", a[0].len(), a.len());
    }
    if a.len() != b.len() {
        panic!("err a.len():{:?} != b.len():{:?}", a.len(), b.len());
    }

    let rows = a.len();
    let cols = a.len();

    let mut c: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); cols];
    for i in 0..rows {
        for j in 0..cols {
            c[i] += a[i][j] * b[j];
        }
    }
    c
}

fn add_vv(a: Vec<Complex64>, b: Vec<Complex64>) -> Vec<Complex64> {
    if a.len() != b.len() {
        panic!("err a.len():{:?} != b.len():{:?}", a.len(), b.len());
    }
    let mut c: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
    c
}

// mul_vv_el multiplies elements of one vector by the elements of another vector
fn mul_vv_el(a: Vec<Complex64>, b: Vec<Complex64>) -> Vec<Complex64> {
    if a.len() != b.len() {
        panic!("err a.len():{:?} != b.len():{:?}", a.len(), b.len());
    }
    let mut c: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] * b[i];
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate rand;
    use rand::Rng;

    #[test]
    fn test_dft_simple_values() {
        let values: Vec<f64> = vec![0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
        let r = dft(&values);
        assert_eq!(r.len(), 8);

        assert_eq!(format!("{:.2}", r[0]), "3.70+0.00i");
        assert_eq!(format!("{:.2}", r[1]), "-0.30-0.97i");
        assert_eq!(format!("{:.2}", r[2]), "-0.30-0.40i");
        assert_eq!(format!("{:.2}", r[3]), "-0.30-0.17i");
        assert_eq!(format!("{:.2}", r[4]), "-0.30+0.00i");

        // expect result similar to initial values
        let o = idft(&r);
        assert_eq!(format!("{:.1}", o[0]), "0.2");
        assert_eq!(format!("{:.1}", o[1]), "0.2");
        assert_eq!(format!("{:.1}", o[2]), "0.3");
        assert_eq!(format!("{:.1}", o[3]), "0.4");
        assert_eq!(format!("{:.1}", o[4]), "0.5");
        assert_eq!(format!("{:.1}", o[5]), "0.6");
        assert_eq!(format!("{:.1}", o[6]), "0.7");
        assert_eq!(format!("{:.1}", o[7]), "0.8");
    }

    #[test]
    fn test_fft_simple_values() {
        let values: Vec<f64> = vec![0.2, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
        let r = fft(&values);

        assert_eq!(r.len(), 8);

        assert_eq!(format!("{:.2}", r[0]), "3.70+0.00i");
        assert_eq!(format!("{:.2}", r[1]), "-0.30-0.97i");
        assert_eq!(format!("{:.2}", r[2]), "-0.30-0.40i");
        assert_eq!(format!("{:.2}", r[3]), "-0.30-0.17i");
        assert_eq!(format!("{:.2}", r[4]), "-0.30+0.00i");

        // expect result similar to initial values
        let o = ifft(&r);
        println!("{:?}", o);
        assert_eq!(format!("{:.1}", o[0]), "0.2");
        assert_eq!(format!("{:.1}", o[1]), "0.2");
        assert_eq!(format!("{:.1}", o[2]), "0.3");
        assert_eq!(format!("{:.1}", o[3]), "0.4");
        assert_eq!(format!("{:.1}", o[4]), "0.5");
        assert_eq!(format!("{:.1}", o[5]), "0.6");
        assert_eq!(format!("{:.1}", o[6]), "0.7");
        assert_eq!(format!("{:.1}", o[7]), "0.8");
    }

    #[test]
    fn test_dft_random_values() {
        let values: Vec<f64> = rand::thread_rng()
            .sample_iter(rand::distributions::Standard)
            .take(1024)
            .collect();
        let r = dft(&values);
        println!("{:?}", r.len());
        let o = idft(&r);
    }

    #[test]
    fn test_fft_random_values() {
        let values: Vec<f64> = rand::thread_rng()
            .sample_iter(rand::distributions::Standard)
            .take(1024)
            .collect();
        let r = fft(&values);
        println!("{:?}", r.len());
        let o = ifft(&r);
    }
}
