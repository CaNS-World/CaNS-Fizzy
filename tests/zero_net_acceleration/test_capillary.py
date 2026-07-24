#!/usr/bin/env python
import numpy as np


def test_capillary_and_gravity_balance():
    balance = np.loadtxt("forcing_balance.out")
    wall = balance[:, 1:4]
    gravity = balance[:, 4:7]
    surface = balance[:, 7:10]
    boussinesq = balance[:, 10:13]
    applied = balance[:, 13:16]
    residual = balance[:, 16:19]

    np.testing.assert_allclose(wall, 0.0, rtol=0.0, atol=1e-20)
    np.testing.assert_allclose(boussinesq, 0.0, rtol=0.0, atol=1e-20)
    assert np.max(np.abs(gravity)) < 1e-7
    assert np.max(np.abs(surface)) > 1e-4
    np.testing.assert_allclose(
        applied, wall - gravity - surface - boussinesq, rtol=2e-14, atol=2e-15
    )
    # The force controller closes its specified source balance.  The measured
    # residual must still expose the non-conservative pressure contribution of
    # the constant-coefficient split at unequal density.
    assert np.max(np.abs(residual)) > 1e-4


if __name__ == "__main__":
    test_capillary_and_gravity_balance()
