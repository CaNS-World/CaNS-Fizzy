#!/usr/bin/env python
import numpy as np


def test_marangoni_balance():
    balance = np.loadtxt("forcing_balance.out")
    wall = balance[:, 1:4]
    gravity = balance[:, 4:7]
    surface = balance[:, 7:10]
    boussinesq = balance[:, 10:13]
    applied = balance[:, 13:16]
    residual = balance[:, 16:19]
    marangoni = balance[:, 19:22]

    # The force of the variable surface tension has a non-zero mean, which the
    # controller balances together with the other forces; at equal density
    # the total mixture momentum then stays constant.
    assert np.max(np.abs(marangoni)) > 1e-4
    np.testing.assert_allclose(
        applied, wall - gravity - surface - boussinesq - marangoni, rtol=2e-13, atol=2e-14
    )
    np.testing.assert_allclose(residual, 0.0, rtol=0.0, atol=1e-11)


if __name__ == "__main__":
    test_marangoni_balance()
