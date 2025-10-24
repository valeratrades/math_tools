use clap::{Parser, Subcommand};
use color_eyre::eyre::Result;
use math_tools::polynomial::solve_quadratic;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
	#[command(subcommand)]
	command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
	/// Solve polynomial equations
	Polynomial {
		/// Coefficient for x^2 term (a*x^2)
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		a: f64,

		/// Coefficient for x^1 term (b*x)
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		b: f64,

		/// Coefficient for x^0 term (constant c)
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		c: f64,

		/// Coefficient for x^3 term (d*x^3) - not yet implemented
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		d: f64,

		/// Coefficient for x^4 term (e*x^4) - not yet implemented
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		e: f64,

		/// Coefficient for x^5 term (f*x^5) - not yet implemented
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		f: f64,

		/// Coefficient for x^6 term (g*x^6) - not yet implemented
		#[arg(short, long, default_value = "0.0", allow_hyphen_values = true)]
		g: f64,
	},
}

fn main() -> Result<()> {
	color_eyre::install()?;

	let cli = Cli::parse();

	match cli.command {
		Commands::Polynomial { a, b, c, d, e, f, g } => {
			// Check which is the highest non-zero coefficient to determine the degree
			let coefficients = vec![(g, 6, "g"), (f, 5, "f"), (e, 4, "e"), (d, 3, "d"), (a, 2, "a"), (b, 1, "b"), (c, 0, "c")];

			// Find the highest degree (first non-zero coefficient)
			let degree = coefficients.iter().find(|(coef, _, _)| coef.abs() > f64::EPSILON).map(|(_, deg, _)| *deg).unwrap_or(0);

			println!("Polynomial equation:");
			print_polynomial(&coefficients);

			match degree {
				0 =>
					if c.abs() < f64::EPSILON {
						println!("Trivial equation: 0 = 0 (infinitely many solutions)");
					} else {
						println!("No solutions: {c} = 0 is false");
					},
				1 => {
					// Linear equation: bx + c = 0
					let x = -c / b;
					println!("Linear equation solution: x = {x}");
				}
				2 => {
					// Quadratic equation: ax^2 + bx + c = 0
					solve_quadratic(a, b, c)?;
				}
				_ => {
					println!("Polynomials of degree {degree} are not yet implemented.");
					println!("Currently only quadratic equations (degree 2) are supported.");
				}
			}
		}
	}

	Ok(())
}

fn print_polynomial(coefficients: &[(f64, usize, &str)]) {
	let mut first = true;
	for (coef, degree, _name) in coefficients.iter() {
		if coef.abs() < f64::EPSILON {
			continue;
		}

		if !first && *coef > 0.0 {
			print!(" + ");
		} else if *coef < 0.0 {
			if first {
				print!("-");
			} else {
				print!(" - ");
			}
		}

		let abs_coef = coef.abs();

		match *degree {
			0 => print!("{abs_coef}"),
			1 =>
				if (abs_coef - 1.0).abs() < f64::EPSILON {
					print!("x");
				} else {
					print!("{abs_coef}x");
				},
			_ =>
				if (abs_coef - 1.0).abs() < f64::EPSILON {
					print!("x^{degree}");
				} else {
					print!("{abs_coef}x^{degree}");
				},
		}

		first = false;
	}

	if first {
		print!("0");
	}

	println!(" = 0");
}
