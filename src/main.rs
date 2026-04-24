//mod grid;
mod integration;
mod interpolation;

//use grid::make_grid;
use integration::*;
use interpolation::*;

// We want to solve
// y(x) = 1 + 0.5\int_{-1}^1 x * t * y(t) dt
//const NZ: usize = 100;
const NUMBER_OF_POLS: usize = 50;
//const INF_BOUND: f64 = 0.001;
//const SUP_BOUND: f64 = 7.0;

fn main() {
    //let x_grid = make_grid(NZ, INF_BOUND, SUP_BOUND); //x_grid.len() = 2 * NZ 

    let (x, w) = gauss_legendre_roots_weights(NUMBER_OF_POLS);
    println!("Finished glrw");
    let f = vec![1.0; NUMBER_OF_POLS];

    let result = quadrature_integrate(&w, &f);

    println!("Result is {result}");

    let barycentric_weights = barycentric_weights(&x);
    let point = 0.37;
    let y = barycentric_interpolate(&x, &f, &barycentric_weights, point);

    println!("The interpolated point is {}", y);
}
