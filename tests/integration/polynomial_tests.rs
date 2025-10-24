use std::process::Command;

/// Test runner that executes polynomial and captures output
fn run_polynomial(args: &[&str]) -> String {
	let output = Command::new("cargo")
		.args(["run", "--quiet", "--"])
		.args(["polynomial"])
		.args(args)
		.env("RUSTC_WRAPPER", "")
		.output()
		.expect("Failed to execute command");

	String::from_utf8_lossy(&output.stdout).to_string()
}

#[test]
fn test_quadratic_exact_integers() {
	let output = run_polynomial(&["-a", "1", "-b", "-5", "-c", "6"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^2 - 5x + 6 = 0

	Two distinct roots:
	  x ∈ {3, 2}
	  ≈ {3, 2}
	");
}

#[test]
fn test_quadratic_with_sqrt() {
	let output = run_polynomial(&["-a", "1", "-b", "-2", "-c", "-1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^2 - 2x - 1 = 0

	Two distinct roots:
	  x ∈ {1 + sqrt(2), 1 - sqrt(2)}
	  ≈ {2.414213562373095, -0.41421356237309515}
	");
}

#[test]
fn test_quadratic_repeated_root() {
	let output = run_polynomial(&["-a", "1", "-b", "2", "-c", "1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^2 + 2x + 1 = 0

	One repeated root: x = -1
	");
}

#[test]
fn test_quadratic_complex_roots() {
	let output = run_polynomial(&["-a", "1", "-b", "1", "-c", "1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^2 + x + 1 = 0

	Two complex conjugate roots:
	  x1 = -0.5 + 0.8660254037844386i
	  x2 = -0.5 - 0.8660254037844386i
	");
}

#[test]
fn test_quadratic_simplified_sqrt() {
	let output = run_polynomial(&["-a", "1", "-b", "4", "-c", "1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^2 + 4x + 1 = 0

	Two distinct roots:
	  x ∈ {-2 + sqrt(3), -2 - sqrt(3)}
	  ≈ {-0.2679491924311228, -3.732050807568877}
	");
}

#[test]
fn test_cubic_perfect_cube() {
	let output = run_polynomial(&["-a", "1", "-d", "-8"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^3 - 8 = 0

	One real root:
	  x = 2
	Two complex conjugate roots:
	  (complex cube roots of 8)
	");
}

#[test]
fn test_cubic_general_cbrt() {
	let output = run_polynomial(&["-a", "1", "-d", "-10"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^3 - 10 = 0

	One real root:
	  x = cbrt(10)
	  x ≈ 2.154434690031884
	Two complex conjugate roots:
	  (complex cube roots of 10)
	");
}

#[test]
fn test_cubic_triple_root() {
	let output = run_polynomial(&["-a", "1", "-b", "-3", "-c", "3", "-d", "-1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^3 - 3x^2 + 3x - 1 = 0

	Triple root: x = 1
	");
}

#[test]
fn test_cubic_three_real_roots() {
	let output = run_polynomial(&["-a", "1", "-b", "-6", "-c", "11", "-d", "-6"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^3 - 6x^2 + 11x - 6 = 0

	Three distinct real roots:
	  x ≈ {3, 2, 0.9999999999999998}
	");
}

#[test]
fn test_quartic_perfect_fourth() {
	let output = run_polynomial(&["-a", "1", "-e", "-16"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^4 - 16 = 0

	Four roots:
	  x ∈ {2, -2}  (real)
	  x ∈ {2i, -2i}  (imaginary)
	");
}

#[test]
fn test_quartic_general_fourth_root() {
	let output = run_polynomial(&["-a", "1", "-e", "-5"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^4 - 5 = 0

	Four roots:
	  x = ±sqrt(sqrt(5))  (real)
	  x = ±i*sqrt(sqrt(5))  (imaginary)
	  x ≈ ±1.4953487812212205, ±1.4953487812212205i
	");
}

#[test]
fn test_quartic_biquadratic() {
	let output = run_polynomial(&["-a", "1", "-e", "-1"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	x^4 - 1 = 0

	Four roots:
	  x ∈ {1, -1}  (real)
	  x ∈ {1i, -1i}  (imaginary)
	");
}

#[test]
fn test_linear() {
	let output = run_polynomial(&["-b", "2", "-c", "-10"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	2x - 10 = 0

	Linear solution: x = 5
	");
}

#[test]
fn test_dynamic_letter_interpretation() {
	let output = run_polynomial(&["-a", "4", "-d", "-10"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	4x^3 - 10 = 0

	One real root:
	  x ≈ 1.3572088082974534
	Two complex conjugate roots:
	  x ≈ -0.6786044041487267 ± 1.1753773062255988i
	");
}

#[test]
fn test_constant_zero() {
	let output = run_polynomial(&["-c", "0"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	0 = 0

	Trivial equation: 0 = 0 (infinitely many solutions)
	");
}

#[test]
fn test_constant_nonzero() {
	let output = run_polynomial(&["-c", "5"]);
	insta::assert_snapshot!(output, @r"
	Polynomial equation:
	5 = 0

	No solutions: 5 = 0 is false
	");
}
