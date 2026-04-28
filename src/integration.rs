const PI: f64 = std::f64::consts::PI;


/// Computes Gauss–Legendre quadrature nodes and weights on the interval [-1, 1].
///
/// This function constructs the Gauss–Legendre integration rule of order n,
/// producing nodes (roots) and weights used for numerical integration.
///
/// The nodes are computed by solving for the roots of the Legendre polynomial
/// \( P_n(x) \) using a Newton–Raphson iteration. Polynomial values and derivatives
/// are evaluated via a stable three-term recurrence relation.
///
/// Symmetry of the Legendre roots is exploited: only half of the roots are computed
/// explicitly, and the remainder are obtained by reflection.
///
/// # Arguments
/// * n - Number of quadrature points (order of the rule)
///
/// # Returns
/// A tuple (x, w):
/// * x - Quadrature nodes in [-1, 1]
/// * w - Corresponding quadrature weights
///
/// # Complexity
/// * O(n**2) per root refinement step (due to recurrence evaluation)
///
/// # Remarks
/// * Suitable for moderate values of n
/// * For large n, eigenvalue-based methods (Golub–Welsch) are more stable and efficient
///
/// # References
/// * Numerical Recipes (Gauss–Legendre routine)
/// * Szegő — Orthogonal Polynomials
pub fn gauss_legendre_roots_weights(
    n: usize,
) -> (Vec<f64>, Vec<f64>){
    
    let tol = 1e-14_f64;

    let x1 = -1.0;
    let x2 = 1.0;

    let m = (n + 1)/2;
    let xm = 0.5 * ( x2 + x1 );
    let xl = 0.5 * ( x2 - x1 );
    let mut z: f64;
    let mut z1: f64;

    let mut p1: f64;
    let mut p2: f64;
    let mut p3: f64;
    let mut pp: f64;

    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    
    for i in 1..=m {
        // Initial guess for the root

        z = f64::cos(PI * ((i as f64) - 0.25) / ( (n as f64) + 0.5) );
        'calculation: loop {
            p1 = 1.0;
            p2 = 0.0;
            for j in 1..=n {
                p3 = p2;
                p2 = p1;
                p1 = ((2.0*(j as f64) - 1.0) * z * p2 - ((j as f64) - 1.0)*p3)/(j as f64);
            }
            //Compute the derivative of the legendre polynomial
            pp = (n as f64)*(z*p1 - p2)/(z*z - 1.0);
            z1 = z;
            z = z1 - p1/pp;

            if (z-z1).abs() < tol{
                break 'calculation;
            }
        }
        
        x[i-1] = xm - xl * z;
        x[n - 1 - (i-1)] = xm + xl * z;

        w[i-1] = 2.0 * xl / ( (1.0-z*z)*pp*pp);
        w[n - 1  - (i-1)] = w[i-1];
    }

    return (x, w);

}

/// Computes the weighted sum used in numerical quadrature.
///
/// This function evaluates the discrete integral approximation:
/// \[ \sum_i w_i \cdot f_i \]
/// where w are the quadrature weights and f are the function
/// values evaluated at the corresponding nodes.
///
/// # Arguments
/// * w - Quadrature weights
/// * f - Function values at the quadrature nodes
///
/// # Returns
/// The approximated integral value.
///
/// # Complexity
/// * O(n)
///
/// # Remarks
/// * Assumes w[i] corresponds to f[i]
/// * No validation is performed beyond slice bounds
pub fn quadrature_integrate(
    w: &[f64],
    f: &[f64],
) -> f64 {
    assert_eq!(w.len(), f.len());
    let mut acc = 0.0;
    for i in 0..w.len() {
        acc += w[i] * f[i];
    }
    acc  
}
