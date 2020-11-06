use num::complex::{Complex, Complex64};
use std::f64::consts::PI;

// dft computes the Discrete Fourier Transform
pub fn dft(x: &Vec<f64>) -> Vec<Complex64> {
    let mut x_compl: Vec<Complex64> = vec![Complex::new(0_f64, 0_f64); x.len()];
    for i in 0..x.len() {
        x_compl[i] = Complex::new(x[i], 0_f64);
    }

    let mut w = Complex::new(0_f64, -2_f64 * PI / x.len() as f64);

    // f_k = SUM{n=0, N-1} f_n * e^(-j2pi*k*n)/N
    // https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    let mut f: Vec<Vec<Complex64>> = Vec::new();
    for i in 0..x.len() {
        let mut f_k: Vec<Complex64> = Vec::new();
        for j in 0..x.len() {
            let i_compl = Complex::new(0_f64, i as f64);
            let j_compl = Complex::new(0_f64, j as f64);
            let fe = (w * i_compl * j_compl).exp();
            f_k.push(fe);
        }
        f.push(f_k.clone());
    }
    let r = mul_vv(f, x_compl);
    r
}

// mul_vv multiplies a Matrix by a Vector
fn mul_vv(a: Vec<Vec<Complex64>>, b: Vec<Complex64>) -> Vec<Complex64> {
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

#[cfg(test)]
mod tests {
    use super::*;

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
    }
}
