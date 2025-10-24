## polynomial
degrees are automatically inferred from largest provided letter, so to solve eg `4x^3 - 10 = 0`, you'd do:
```sh
math_tools polynomial -a 4 -d "-10"
```
outputting:
```
Polynomial equation:
4x^3 - 10 = 0

One real root:
  x ≈ 1.3572088082974534
Two complex conjugate roots:
  x ≈ -0.6786044041487267 ± 1.1753773062255988i
```

For more examples, look at [integration tests](../tests/integration/polynomial_tests.rs)
