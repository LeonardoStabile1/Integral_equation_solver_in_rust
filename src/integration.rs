use gauss_quad::GaussLegendre;

// u(x) = 1 + 0.5 * \int_{-1}^{1} x * t * u(t) dt

pub fn integrand(
    x: f64,
    old_f: &Vec<f64>,
    number_of_pols: usize,
) -> f64 {
    let gl = GaussLegendre::init(number_of_pols);

    x * gl.integrate(-1.0, 1.0, |t: f64| {
        let idx = gl.nodes.iter()
            .position(|&node| node == t)
            .unwrap();
        t * old_f[idx]
    })
}