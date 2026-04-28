// Compute barycentric weights
pub fn barycentric_weights(
    x: &[f64],
) -> Vec<f64> {
    let n = x.len();
    let mut w = vec![1.0; n];

    for j in 0..n {
        for k in 0..n{
            if j!= k {
                w[j] /= x[j] - x[k];
            }
        }
    }

    w
}

// Interpolation method
pub fn barycentric_interpolate(
    x_nodes: &[f64],
    f: &[f64],
    w: &[f64],
    x: f64
) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for j in 0..x_nodes.len() {
        let dx = x - x_nodes[j];

        if dx.abs() < 1e-14{
            return f[j];
        }
        num += w[j]*f[j] / dx;
        den += w[j] / dx;
    }

    return num / den;
}
