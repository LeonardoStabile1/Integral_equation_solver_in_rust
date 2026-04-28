//mod grid;
mod integration;
//mod interpolation;
use integration::*;
use nalgebra::DVector;



// We want to solve
// y(x) = 1 + 0.5\int_{-1}^1 y(t) dt
const NUMBER_OF_POLS: usize = 50;


fn main() {

    let (x, w) = gauss_legendre_roots_weights(NUMBER_OF_POLS);
    let f = vec![1.0; NUMBER_OF_POLS];
    let mut g = vec![0.0; NUMBER_OF_POLS];

    g = solve(&f, &w, &x);

    println!("{:?}", g)

}

fn solve(f: &[f64], w: &[f64], x: &[f64]) -> Vec<f64>{
    let mut next_f = vec![0.0; f.len()];
    let mut old_f = f.to_vec();
    let tolerance = 1e-5_f64;
    'solver: loop{
        next_f = integrand(&old_f, &w, &x);
        let error = relative_error(&old_f, &next_f);
        if error < tolerance{
            break 'solver;
        }
        println!("Error: {error}");
        std::mem::swap(&mut old_f, &mut next_f);
    }

    next_f
}

fn relative_error(f: &[f64], new_f: &[f64]) -> f64{
    assert_eq!(f.len(), new_f.len());
    let vf = DVector::from_column_slice(f);
    let vnew = DVector::from_column_slice(new_f);

    let num = (&vf - &vnew).norm();
    let den = vf.norm();

    num / den
}

fn integrand(f: &[f64], w: &[f64], x: &[f64]) -> Vec<f64>{
    let mut integrand = vec![0.0; NUMBER_OF_POLS];
    let mut new_f = vec![0.0; NUMBER_OF_POLS];
    let mut integral_value: f64;

    for i in 0..NUMBER_OF_POLS{
        for j in 0..NUMBER_OF_POLS {
            integrand[j] = (x[j].powi(2) + x[i].powi(2))*f[j];
        }
        integral_value = quadrature_integrate(&w, &integrand);

        new_f[i] = 1.0 + 0.1*integral_value;
    }
    new_f
}
