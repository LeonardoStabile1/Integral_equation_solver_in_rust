pub fn make_grid(
    nx: usize,
    inf_bound: f64,
    sup_bound: f64,
) -> Vec<f64> {
    assert!(inf_bound > 0.0);
    assert!(nx > 1);

    println!("Building a log-scale grid of {nx} poins between [{inf_bound}, {sup_bound}]");

    let log_inf = inf_bound.log10();
    let log_sup = sup_bound.log10();

    let step = (log_sup - log_inf) / (nx as f64 - 1.0);

    let mut xvec = Vec::with_capacity(2*nx);

    for i in 0..nx {
        let value = log_inf + step * i as f64;
        xvec.push(10.0_f64.powf(value));
        xvec.push(-10.0_f64.powf(value));
    }

    xvec
}
