mod grid;
mod integration;

use grid::make_grid;
use integration::*;

// We want to solve
// y(x) = 1 + 0.5\int_{-1}^1 x * t * u(t) dt

fn main() {

    const NZ: usize = 100;
    const NUMBER_OF_POLS: usize = 100;

    let x_grid = make_grid(NZ, 0.001, 7.0); //x_grid.len() = 2 * NZ 

    let (x, w) = gauss_legendre_roots_weights(10);

    let f = vec![1.0; 10];

    let result = quadrature_integrate(&w, &f);

    println!("Result is {result}");
}
