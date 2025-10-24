use std::collections::BTreeMap;

use clap::{Parser, Subcommand};
use color_eyre::eyre::Result;
use math_tools::polynomial::solve_polynomial;

#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
	#[command(subcommand)]
	command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
	/// Solve polynomial equations
	Polynomial {
		/// Coefficient a (degree determined by highest non-zero letter)
		#[arg(short, long, allow_hyphen_values = true)]
		a: Option<f64>,

		/// Coefficient b
		#[arg(short, long, allow_hyphen_values = true)]
		b: Option<f64>,

		/// Coefficient c
		#[arg(short, long, allow_hyphen_values = true)]
		c: Option<f64>,

		/// Coefficient d
		#[arg(short, long, allow_hyphen_values = true)]
		d: Option<f64>,

		/// Coefficient e
		#[arg(short, long, allow_hyphen_values = true)]
		e: Option<f64>,

		/// Coefficient f
		#[arg(short, long, allow_hyphen_values = true)]
		f: Option<f64>,

		/// Coefficient g
		#[arg(short, long, allow_hyphen_values = true)]
		g: Option<f64>,
	},
}

fn main() -> Result<()> {
	color_eyre::install()?;

	let cli = Cli::parse();

	match cli.command {
		Commands::Polynomial { a, b, c, d, e, f, g } => {
			// Build a BTreeMap of provided coefficients
			let mut coeffs_map = BTreeMap::new();
			if let Some(val) = a {
				coeffs_map.insert('a', val);
			}
			if let Some(val) = b {
				coeffs_map.insert('b', val);
			}
			if let Some(val) = c {
				coeffs_map.insert('c', val);
			}
			if let Some(val) = d {
				coeffs_map.insert('d', val);
			}
			if let Some(val) = e {
				coeffs_map.insert('e', val);
			}
			if let Some(val) = f {
				coeffs_map.insert('f', val);
			}
			if let Some(val) = g {
				coeffs_map.insert('g', val);
			}

			if coeffs_map.is_empty() {
				println!("No coefficients provided. Please provide at least one coefficient.");
				return Ok(());
			}

			// Find the highest letter provided to determine n
			let highest_letter = *coeffs_map.keys().last().unwrap();
			let n = (highest_letter as u8 - b'a') as usize;

			// Build coefficient vector where coeffs[i] is the coefficient of x^(n-i)
			// So for `-a 4 -d -10`, we have highest_letter='d' (n=3)
			// 'a' maps to x^3, 'b' maps to x^2, 'c' maps to x^1, 'd' maps to x^0
			let mut coeffs = vec![0.0; n + 1];
			for (letter, value) in coeffs_map.iter() {
				let letter_index = (*letter as u8 - b'a') as usize;
				coeffs[letter_index] = *value;
			}

			// Print the polynomial equation
			print_polynomial(&coeffs, n);

			// Solve the polynomial
			solve_polynomial(&coeffs)?;
		}
	}

	Ok(())
}

fn print_polynomial(coeffs: &[f64], n: usize) {
	println!("Polynomial equation:");
	let mut first = true;

	for (i, &coef) in coeffs.iter().enumerate() {
		if coef.abs() < f64::EPSILON {
			continue;
		}

		let degree = n - i;

		if !first && coef > 0.0 {
			print!(" + ");
		} else if coef < 0.0 {
			if first {
				print!("-");
			} else {
				print!(" - ");
			}
		}

		let abs_coef = coef.abs();

		match degree {
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

	println!(" = 0\n");
}
