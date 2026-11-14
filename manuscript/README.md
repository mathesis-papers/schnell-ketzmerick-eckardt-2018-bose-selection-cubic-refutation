# Two cubic baths and an exponent-dependent selection bound

The note addresses the non-Ohmic expectation in Schnell, Ketzmerick and Eckardt, arXiv:1802.04088v2, Section V, page 11, equation (55) and its surrounding paragraph. Two cubic product-coupled baths have five uniquely selected modes on an open set. The same five selected modes persist after adding any finite number of sufficiently weakly coupled higher levels.

For one or more baths with a common positive odd integer exponent `q`, positive coupling weights and nonzero signed strengths, the note proves `M_S <= (q+1) M_B - 1` when the selected principal matrix has a one-dimensional kernel. This extends the Ohmic bound of Schnell, Ketzmerick and Eckardt (2018). For two cubic baths the upper bound is seven under that hypothesis; attaining seven is not established by these arguments.

## Files and reproduction

The source archive contains `note.tex`, this `README.md`, `metadata.json`, `verify.py`, `verification.json` and `MANIFEST.sha256`. The typeset `note.pdf` accompanies it separately. The bibliography is inside `note.tex`; external data and BibTeX are unnecessary.

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
python3 verify.py > verification.json
sha256sum -c MANIFEST.sha256
```

Alternatively, use a standard LaTeX distribution:

```sh
pdflatex -no-shell-escape note.tex
pdflatex -no-shell-escape note.tex
python3 verify.py > verification.json
```

The Python program uses only the standard library and requires Python 3.9 or later. Its expected JSON ends with `"checks": "pass"`; the output SHA256 is `626653c1609c72a0c824a9205f376d461ed7193c7f45e5d87911dfacfd919edf`.

## Mathematical coverage

At energies `(1,2,3,4,5)`, the normalized coupling weights are `(5,4,4,4,4)/21` and `(6,6,3,3,3)/21`, and the occupation fractions are `(14,56,21,8,1)/100`. The exact arithmetic reproduces the integer-scaled matrix, zero kernel residual, signed Pfaffians `(504,2016,756,288,36)`, rank four and principal determinant `254016`. Five polynomial Pfaffian-kernel identities are checked by coefficient cancellation. Additional rational checks give exactly five selected modes in dimensions six, seven and eight, and verify the binomial matrix factorizations at exponents one, three and five. These finite examples illustrate, rather than replace, the general written proofs.

The declaration `D5.S3.StatisticalMechanics.BoseSelectionCubicTwoBathRefutation.result` in the [pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/6ffa76e4424ce4c0d1048e56afb0724762ed8285/D5/S3/StatisticalMechanics/BoseSelectionCubicTwoBathRefutation.lean) negates the all-exponent, all-dimension empty-interior claim. Its construction gives a nonempty open violating set at exponent three in dimension five, using full selection. The PDF also links the pinned normalization and signed-power definitions. Uniqueness, positive transition rates, exponent stability, higher-level extensions, the odd-monomial upper bound and same-sign-bath selection have written-proof coverage in the note.

The coverage table distinguishes the non-Ohmic counterexample and replacement statement from the source's Ohmic result, the repeated thesis expectation, the standard Pfaffian method and unrelated continuum, random-rate, driven-system and interaction questions. Knebel et al. (2015) supplies the Pfaffian method; Geiger et al. (2018) treats related antisymmetric networks; Petiziol and Eckardt (2024) gives a different five-mode construction using frequency-selective reservoirs. No priority assertion is made.

Haobo Ma (ChronoAI Pte Ltd, Singapore); Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore).

Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
