mod grid;
mod integration;

use grid::make_grid;
use integration::integrand;

// We want to solve
// y(x) = 1 + 0.5\int_{-1}^1 x * t * u(t) dt

fn main() {

    const NZ: usize = 100;
    const NUMBER_OF_POLS: usize = 100;

    let x_grid = make_grid(NZ, 0.001, 7.0); //x_grid.len() = 2 * NZ 

    let mut x_val: f64;

    let function_values: Vec<f64> = vec![0.0; 2*NZ];

    let mut result: f64;

    for i in 0..x_grid.len(){
        x_val = x_grid[i];
        result = integrand(x_val, &function_values, NUMBER_OF_POLS);
        println!("{result}");
    }
}