use color_eyre::eyre::{Result, bail};

/// Helper to simplify a fraction by finding GCD
fn gcd(mut a: i64, mut b: i64) -> i64 {
	a = a.abs();
	b = b.abs();
	while b != 0 {
		let t = b;
		b = a % b;
		a = t;
	}
	a
}

/// Helper to check if a number is close to an integer
fn is_near_integer(x: f64) -> Option<i64> {
	let rounded = x.round();
	if (x - rounded).abs() < 1e-9 { Some(rounded as i64) } else { None }
}

/// Helper to simplify sqrt(n) = a*sqrt(b) where b is square-free
fn simplify_sqrt(n: f64) -> (i64, i64) {
	if n < 0.0 {
		return (0, 0);
	}

	let n_int = n.round() as i64;
	if (n - n_int as f64).abs() > 1e-9 {
		// Not an integer, can't simplify nicely
		return (1, n_int);
	}

	let mut outside = 1i64;
	let mut inside = n_int;

	let mut i = 2i64;
	while i * i <= inside {
		while inside % (i * i) == 0 {
			outside *= i;
			inside /= i * i;
		}
		i += 1;
	}

	(outside, inside)
}

/// Format an exact quadratic solution
fn format_exact_quadratic(b: f64, discriminant: f64, a: f64, plus: bool) -> String {
	// Check if coefficients are near integers
	let b_int = is_near_integer(b);
	let a_int = is_near_integer(a);
	let disc_int = is_near_integer(discriminant);

	if let (Some(b), Some(a), Some(disc)) = (b_int, a_int, disc_int) {
		if disc < 0 {
			return String::new(); // Complex, handled elsewhere
		}

		let (sqrt_coef, sqrt_inside) = simplify_sqrt(disc as f64);

		// Check if it's a perfect square
		if sqrt_inside == 1 {
			let numerator = if plus { -b + sqrt_coef } else { -b - sqrt_coef };
			let denominator = 2 * a;
			let g = gcd(numerator, denominator);
			let num = numerator / g;
			let den = denominator / g;

			if den == 1 {
				return format!("{}", num);
			} else {
				return format!("{}/{}", num, den);
			}
		}

		// Format with square root
		let sign = if plus { "+" } else { "-" };
		let numerator_const = -b;
		let denominator = 2 * a;

		// Try to simplify the fraction
		if numerator_const == 0 {
			// Just the square root part
			let g = gcd(sqrt_coef, denominator);
			let sqrt_num = sqrt_coef / g;
			let den = denominator / g;

			if sqrt_num == 1 {
				if den == 1 {
					return format!("{}sqrt({})", if plus { "" } else { "-" }, sqrt_inside);
				} else {
					return format!("{}sqrt({})/{}", if plus { "" } else { "-" }, sqrt_inside, den);
				}
			} else {
				if den == 1 {
					return format!("{}{}*sqrt({})", if plus { "" } else { "-" }, sqrt_num, sqrt_inside);
				} else {
					return format!("{}{}*sqrt({})/{}", if plus { "" } else { "-" }, sqrt_num, sqrt_inside, den);
				}
			}
		} else {
			// Both constant and square root part
			let g = gcd(gcd(numerator_const.abs(), sqrt_coef), denominator);
			let const_num = numerator_const / g;
			let sqrt_num = sqrt_coef / g;
			let den = denominator / g;

			let sqrt_str = if sqrt_num == 1 {
				format!("sqrt({})", sqrt_inside)
			} else {
				format!("{}*sqrt({})", sqrt_num, sqrt_inside)
			};

			let result = if den == 1 {
				format!("{} {} {}", const_num, sign, sqrt_str)
			} else {
				format!("({} {} {})/{}", const_num, sign, sqrt_str, den)
			};

			return result;
		}
	}

	String::new()
}

/// Main entry point for solving polynomial equations of any degree
pub fn solve_polynomial(coeffs: &[f64]) -> Result<()> {
	let degree = coeffs.len() - 1;

	// Find the actual degree (highest non-zero coefficient)
	let mut actual_degree = degree;
	for (i, &c) in coeffs.iter().enumerate() {
		if c.abs() > f64::EPSILON {
			actual_degree = degree - i;
			break;
		}
	}

	match actual_degree {
		0 => solve_constant(coeffs[degree]),
		1 => solve_linear(coeffs[degree - 1], coeffs[degree]),
		2 => solve_quadratic(coeffs[degree - 2], coeffs[degree - 1], coeffs[degree]),
		3 => solve_cubic(coeffs[degree - 3], coeffs[degree - 2], coeffs[degree - 1], coeffs[degree]),
		4 => solve_quartic(coeffs[degree - 4], coeffs[degree - 3], coeffs[degree - 2], coeffs[degree - 1], coeffs[degree]),
		5..=8 => solve_numerical(coeffs, actual_degree),
		_ => {
			bail!("Polynomials of degree {} are not supported", actual_degree);
		}
	}
}

fn solve_constant(c: f64) -> Result<()> {
	if c.abs() < f64::EPSILON {
		println!("Trivial equation: 0 = 0 (infinitely many solutions)");
	} else {
		println!("No solutions: {c} = 0 is false");
	}
	Ok(())
}

fn solve_linear(b: f64, c: f64) -> Result<()> {
	if b.abs() < f64::EPSILON {
		return solve_constant(c);
	}
	let x = -c / b;
	println!("Linear solution: x = {x}");
	Ok(())
}

/// Solve quadratic equation: ax^2 + bx + c = 0
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Result<()> {
	if a.abs() < f64::EPSILON {
		return solve_linear(b, c);
	}

	let discriminant = b * b - 4.0 * a * c;

	if discriminant < -f64::EPSILON {
		// Complex roots
		let real = -b / (2.0 * a);
		let imag = (-discriminant).sqrt() / (2.0 * a);
		println!("Two complex conjugate roots:");
		println!("  x1 = {real} + {imag}i");
		println!("  x2 = {real} - {imag}i");
	} else if discriminant.abs() < f64::EPSILON {
		// One repeated root
		let x = -b / (2.0 * a);
		println!("One repeated root: x = {x}");
	} else {
		// Two distinct real roots
		let sqrt_discriminant = discriminant.sqrt();
		let x1 = (-b + sqrt_discriminant) / (2.0 * a);
		let x2 = (-b - sqrt_discriminant) / (2.0 * a);

		// Try to format exact solutions
		let exact1 = format_exact_quadratic(b, discriminant, a, true);
		let exact2 = format_exact_quadratic(b, discriminant, a, false);

		if !exact1.is_empty() && !exact2.is_empty() {
			println!("Two distinct roots:");
			println!("  x ∈ {{{}, {}}}", exact1, exact2);
			println!("  ≈ {{{}, {}}}", x1, x2);
		} else {
			println!("Two distinct roots:");
			println!("  x1 = {x1}");
			println!("  x2 = {x2}");
		}
	}
	Ok(())
}

/// Solve cubic equation: ax^3 + bx^2 + cx + d = 0 using Cardano's method
fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Result<()> {
	if a.abs() < f64::EPSILON {
		return solve_quadratic(b, c, d);
	}

	// Normalize to x^3 + px + q = 0 form using depressed cubic
	let p = (3.0 * a * c - b * b) / (3.0 * a * a);
	let q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);

	let discriminant = -(4.0 * p * p * p + 27.0 * q * q);

	let offset = -b / (3.0 * a);

	if discriminant > f64::EPSILON {
		// Three distinct real roots (trigonometric method)
		let m = 2.0 * (-p / 3.0).sqrt();
		let theta = (3.0 * q / (p * m)).acos() / 3.0;

		let x1 = m * theta.cos() + offset;
		let x2 = m * (theta - 2.0 * std::f64::consts::PI / 3.0).cos() + offset;
		let x3 = m * (theta - 4.0 * std::f64::consts::PI / 3.0).cos() + offset;

		println!("Three distinct real roots:");
		println!("  x ≈ {{{x1}, {x2}, {x3}}}");
	} else if discriminant.abs() < f64::EPSILON {
		// Multiple roots
		if p.abs() < f64::EPSILON {
			// Triple root at offset
			let offset_int = is_near_integer(offset);
			if let Some(val) = offset_int {
				println!("Triple root: x = {val}");
			} else {
				println!("Triple root: x = {offset}");
			}
		} else {
			let x1 = 3.0 * q / p + offset;
			let x2 = -3.0 * q / (2.0 * p) + offset;

			println!("Multiple roots:");
			println!("  x1 = {x1} (simple)");
			println!("  x2 = {x2} (double)");
		}
	} else {
		// One real root and two complex conjugate roots
		// Check for special case: x^3 = n (pure cube root)
		if b.abs() < f64::EPSILON && c.abs() < f64::EPSILON {
			let n = -d / a;
			let n_int = is_near_integer(n);

			if let Some(n_val) = n_int {
				// Check if it's a perfect cube
				let cbrt_n = (n_val.abs() as f64).cbrt().round() as i64;
				if cbrt_n * cbrt_n * cbrt_n == n_val.abs() {
					let sign = if n_val < 0 { -1 } else { 1 };
					let root = sign * cbrt_n;
					println!("One real root:");
					println!("  x = {root}");
					println!("Two complex conjugate roots:");
					println!("  (complex cube roots of {n_val})");
					return Ok(());
				} else {
					// Not a perfect cube but integer argument
					println!("One real root:");
					if n_val >= 0 {
						println!("  x = cbrt({n_val})");
					} else {
						println!("  x = -cbrt({})", -n_val);
					}
					println!("  x ≈ {}", n.cbrt());
					println!("Two complex conjugate roots:");
					println!("  (complex cube roots of {n_val})");
					return Ok(());
				}
			}
		}

		let sqrt_disc = (-discriminant).sqrt();
		let u = ((-q + sqrt_disc / (3.0 * 3_f64.sqrt())) / 2.0).cbrt();
		let v = ((-q - sqrt_disc / (3.0 * 3_f64.sqrt())) / 2.0).cbrt();

		let x1 = u + v + offset;

		let real = -(u + v) / 2.0 + offset;
		let imag = (u - v) * 3_f64.sqrt() / 2.0;

		println!("One real root:");
		println!("  x ≈ {x1}");
		println!("Two complex conjugate roots:");
		println!("  x ≈ {real} ± {imag}i");
	}

	Ok(())
}

/// Solve quartic equation: ax^4 + bx^3 + cx^2 + dx + e = 0 using Ferrari's method
fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Result<()> {
	if a.abs() < f64::EPSILON {
		return solve_cubic(b, c, d, e);
	}

	// Normalize coefficients
	let b = b / a;
	let c = c / a;
	let d = d / a;
	let e = e / a;

	// Reduce to depressed quartic: y^4 + py^2 + qy + r = 0
	let p = c - 3.0 * b * b / 8.0;
	let q = b * b * b / 8.0 - b * c / 2.0 + d;
	let r = -3.0 * b * b * b * b / 256.0 + c * b * b / 16.0 - b * d / 4.0 + e;

	let offset = -b / 4.0;

	if q.abs() < f64::EPSILON {
		// Biquadratic case: y^4 + py^2 + r = 0
		// Check for special case: x^4 = n (pure fourth root)
		if b.abs() < f64::EPSILON && c.abs() < f64::EPSILON && d.abs() < f64::EPSILON {
			let n = -e / a;
			let n_int = is_near_integer(n);

			if let Some(n_val) = n_int {
				if n_val > 0 {
					// Check if fourth root is an integer
					let fourth_root = (n_val as f64).sqrt().sqrt().round() as i64;
					if fourth_root * fourth_root * fourth_root * fourth_root == n_val {
						println!("Four roots:");
						println!("  x ∈ {{{fourth_root}, -{fourth_root}}}  (real)");
						println!("  x ∈ {{{fourth_root}i, -{fourth_root}i}}  (imaginary)");
						return Ok(());
					} else {
						println!("Four roots:");
						println!("  x = ±sqrt(sqrt({n_val}))  (real)");
						println!("  x = ±i*sqrt(sqrt({n_val}))  (imaginary)");
						println!("  x ≈ ±{}, ±{}i", n.sqrt().sqrt(), n.sqrt().sqrt());
						return Ok(());
					}
				} else if n_val < 0 {
					println!("Four complex roots:");
					println!("  (fourth roots of {n_val})");
					return Ok(());
				}
			}
		}

		// Substitute z = y^2
		let discriminant = p * p - 4.0 * r;
		if discriminant < 0.0 {
			println!("Four complex roots (biquadratic case with complex z)");
			// Would need to implement complex arithmetic for full solution
			return Ok(());
		}

		let z1 = (-p + discriminant.sqrt()) / 2.0;
		let z2 = (-p - discriminant.sqrt()) / 2.0;

		let mut roots = Vec::new();
		let mut exact_strs = Vec::new();

		if z1 > 0.0 {
			let sqrt_z1 = z1.sqrt();
			roots.push(sqrt_z1 + offset);
			roots.push(-sqrt_z1 + offset);

			// Try to format exactly if z1 is a perfect square
			if let Some(z1_int) = is_near_integer(z1) {
				let sqrt_z1_int = (z1_int as f64).sqrt().round() as i64;
				if sqrt_z1_int * sqrt_z1_int == z1_int {
					exact_strs.push(format!("{}", sqrt_z1_int));
					exact_strs.push(format!("-{}", sqrt_z1_int));
				}
			}
		}
		if z2 > 0.0 {
			let sqrt_z2 = z2.sqrt();
			roots.push(sqrt_z2 + offset);
			roots.push(-sqrt_z2 + offset);
		}

		if roots.is_empty() {
			println!("Four complex roots");
		} else {
			if !exact_strs.is_empty() && exact_strs.len() == roots.len() {
				println!("Real roots:");
				print!("  x ∈ {{");
				for (i, s) in exact_strs.iter().enumerate() {
					if i > 0 {
						print!(", ");
					}
					print!("{s}");
				}
				println!("}}");
			} else {
				println!("Real roots:");
				print!("  x ≈ {{");
				for (i, root) in roots.iter().enumerate() {
					if i > 0 {
						print!(", ");
					}
					print!("{root}");
				}
				println!("}}");
			}
		}

		return Ok(());
	}

	// Solve resolvent cubic: 8m^3 + 8pm^2 + (2p^2 - 8r)m - q^2 = 0
	let resolvent_a = 8.0;
	let resolvent_b = 8.0 * p;
	let resolvent_c = 2.0 * p * p - 8.0 * r;
	let resolvent_d = -q * q;

	// Find one real root m of the resolvent cubic (using a simplified approach)
	// For practical purposes, we'll use a numerical method here
	let m = find_cubic_root(resolvent_a, resolvent_b, resolvent_c, resolvent_d);

	// Now solve two quadratics
	let sqrt_2m = (2.0 * m).sqrt();
	let a1 = 1.0;
	let b1 = sqrt_2m;
	let c1 = p / 2.0 + m - q / (2.0 * sqrt_2m);

	let a2 = 1.0;
	let b2 = -sqrt_2m;
	let c2 = p / 2.0 + m + q / (2.0 * sqrt_2m);

	println!("Quartic roots:");
	let mut root_count = 1;

	// Solve first quadratic
	let disc1 = b1 * b1 - 4.0 * a1 * c1;
	if disc1 >= 0.0 {
		let sqrt_disc1 = disc1.sqrt();
		println!("  x{} = {}", root_count, (-b1 + sqrt_disc1) / (2.0 * a1) + offset);
		root_count += 1;
		println!("  x{} = {}", root_count, (-b1 - sqrt_disc1) / (2.0 * a1) + offset);
		root_count += 1;
	}

	// Solve second quadratic
	let disc2 = b2 * b2 - 4.0 * a2 * c2;
	if disc2 >= 0.0 {
		let sqrt_disc2 = disc2.sqrt();
		println!("  x{} = {}", root_count, (-b2 + sqrt_disc2) / (2.0 * a2) + offset);
		root_count += 1;
		println!("  x{} = {}", root_count, (-b2 - sqrt_disc2) / (2.0 * a2) + offset);
	}

	Ok(())
}

/// Find one real root of a cubic equation using Newton-Raphson method
fn find_cubic_root(a: f64, b: f64, c: f64, d: f64) -> f64 {
	// Start with an initial guess
	let mut x = 1.0;

	// Newton-Raphson iteration
	for _ in 0..100 {
		let f = a * x * x * x + b * x * x + c * x + d;
		let f_prime = 3.0 * a * x * x + 2.0 * b * x + c;

		if f_prime.abs() < f64::EPSILON {
			break;
		}

		let x_new = x - f / f_prime;

		if (x_new - x).abs() < 1e-10 {
			return x_new;
		}

		x = x_new;
	}

	x
}

/// Solve polynomial of degree 5-8 using numerical methods (Durand-Kerner)
fn solve_numerical(coeffs: &[f64], degree: usize) -> Result<()> {
	println!("Using numerical method (Durand-Kerner) for degree {degree} polynomial...");

	// Normalize so leading coefficient is 1
	let mut normalized = coeffs.to_vec();
	let leading = coeffs[coeffs.len() - degree - 1];
	for c in normalized.iter_mut() {
		*c /= leading;
	}

	// Initial guesses: points on unit circle
	let mut roots: Vec<(f64, f64)> = Vec::new();
	for k in 0..degree {
		let angle = 2.0 * std::f64::consts::PI * (k as f64) / (degree as f64);
		roots.push((angle.cos() * 0.4 + 0.4, angle.sin() * 0.4));
	}

	// Durand-Kerner iteration
	for iteration in 0..1000 {
		let mut max_change: f64 = 0.0;

		for i in 0..degree {
			// Evaluate polynomial at current root approximation
			let (x, y) = roots[i];
			let (p_val_re, p_val_im) = eval_poly_complex(&normalized, degree, x, y);

			// Compute denominator (product of differences with other roots)
			let mut denom_re = 1.0;
			let mut denom_im = 0.0;

			#[allow(clippy::needless_range_loop)] // we actually do need both `i` and `j` here
			for j in 0..degree {
				if i != j {
					let (xj, yj) = roots[j];
					let diff_re = x - xj;
					let diff_im = y - yj;

					// Multiply complex numbers: denom *= (diff_re + diff_im*i)
					let new_re = denom_re * diff_re - denom_im * diff_im;
					let new_im = denom_re * diff_im + denom_im * diff_re;
					denom_re = new_re;
					denom_im = new_im;
				}
			}

			let denom_mag_sq = denom_re * denom_re + denom_im * denom_im;
			if denom_mag_sq < f64::EPSILON {
				continue;
			}

			let quot_re = (p_val_re * denom_re + p_val_im * denom_im) / denom_mag_sq;
			let quot_im = (p_val_im * denom_re - p_val_re * denom_im) / denom_mag_sq;

			// Update root
			let new_x = x - quot_re;
			let new_y = y - quot_im;

			let change = ((new_x - x) * (new_x - x) + (new_y - y) * (new_y - y)).sqrt();
			max_change = max_change.max(change);

			roots[i] = (new_x, new_y);
		}

		if max_change < 1e-10 {
			println!("Converged after {iteration} iterations");
			break;
		}
	}

	// Print roots
	println!("Numerical roots:");
	for (i, &(re, im)) in roots.iter().enumerate() {
		if im.abs() < 1e-8 {
			println!("  x{} = {}", i + 1, re);
		} else if im > 0.0 {
			println!("  x{} = {} + {}i", i + 1, re, im);
		} else {
			println!("  x{} = {} - {}i", i + 1, re, -im);
		}
	}

	Ok(())
}

/// Evaluate polynomial with complex argument
fn eval_poly_complex(coeffs: &[f64], degree: usize, x: f64, y: f64) -> (f64, f64) {
	let start_idx = coeffs.len() - degree - 1;

	let mut result_re = 0.0;
	let mut result_im = 0.0;

	for (i, &coef) in coeffs[start_idx..].iter().enumerate() {
		let power = degree - i;

		// Compute (x + yi)^power
		let (pow_re, pow_im) = complex_power(x, y, power);

		// Multiply by coefficient and add
		result_re += coef * pow_re;
		result_im += coef * pow_im;
	}

	(result_re, result_im)
}

/// Compute (x + yi)^n using repeated multiplication
fn complex_power(x: f64, y: f64, n: usize) -> (f64, f64) {
	if n == 0 {
		return (1.0, 0.0);
	}

	let mut result_re = 1.0;
	let mut result_im = 0.0;

	for _ in 0..n {
		let new_re = result_re * x - result_im * y;
		let new_im = result_re * y + result_im * x;
		result_re = new_re;
		result_im = new_im;
	}

	(result_re, result_im)
}
