#!/usr/bin/env python
import numpy as np


def test_capillary_momentum_balance():
    balance = np.loadtxt("forcing_balance.out")

    assert np.max(np.abs(balance[:, 7:10])) > 1e-4
    np.testing.assert_allclose(balance[:, 16:19], 0.0, rtol=0.0, atol=1e-11)


if __name__ == "__main__":
    test_capillary_momentum_balance()
