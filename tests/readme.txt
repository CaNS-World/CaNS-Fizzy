this is the directory where tests that are run on GitHub, on push are stored

each test is run from its own directory with `bash testit.sh`, and uses the executable `run/cans` built beforehand
(except `plic_unit`, which builds its own small program, and `single_precision`, which also builds the
code in single precision, in a copy of the source tree, and compares a few cases with double precision);
`testall.sh` builds the variants of the code that the tests need (e.g., ACDI, THINC/QQ) and runs the tests of each one
