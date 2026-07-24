#!/usr/bin/env python
import numpy as np


def test_boussinesq_balance():
    forcing = np.loadtxt("forcing.out")
    balance = np.loadtxt("forcing_balance.out")

    np.testing.assert_allclose(balance[:, 10], -1.0, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(balance[:, 13], 1.0, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(balance[:, 16:19], 0.0, rtol=0.0, atol=2e-12)
    np.testing.assert_allclose(forcing[:, 4:7], 0.0, rtol=0.0, atol=2e-15)


if __name__ == "__main__":
    test_boussinesq_balance()
