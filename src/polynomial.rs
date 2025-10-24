use color_eyre::eyre::{Result, bail};

/// Solve a quadratic equation of the form ax^2 + bx + c = 0
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Result<(Option<f64>, Option<f64>)> {
	// Check if 'a' is zero (not actually a quadratic equation)
	if a.abs() < f64::EPSILON {
		bail!("Coefficient 'a' cannot be zero for a quadratic equation");
	}

	// Calculate discriminant: b^2 - 4ac
	let discriminant = b * b - 4.0 * a * c;

	if discriminant < 0.0 {
		// No real solutions
		println!("No real solutions (discriminant = {discriminant})");
		Ok((None, None))
	} else if discriminant.abs() < f64::EPSILON {
		// One solution (repeated root)
		let x = -b / (2.0 * a);
		println!("One repeated root: x = {x}");
		Ok((Some(x), None))
	} else {
		// Two distinct solutions
		let sqrt_discriminant = discriminant.sqrt();
		let x1 = (-b + sqrt_discriminant) / (2.0 * a);
		let x2 = (-b - sqrt_discriminant) / (2.0 * a);
		println!("Two distinct roots:");
		println!("  x1 = {x1}");
		println!("  x2 = {x2}");
		Ok((Some(x1), Some(x2)))
	}
}
