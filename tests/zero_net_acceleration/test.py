#!/usr/bin/env python
import numpy as np


def test_zero_net_acceleration():
    forcing = np.loadtxt("forcing.out")
    balance = np.loadtxt("forcing_balance.out")

    # The initial midpoint-sampled parabola is normalized to unit bulk velocity.
    nz = 16
    dz = 1.0 / nz
    expected_wall_force = 0.01 * (12.0 - 6.0 * dz) / (1.0 + 0.5 * dz**2)

    np.testing.assert_allclose(-forcing[0, 1], expected_wall_force, rtol=2e-13, atol=2e-15)
    np.testing.assert_allclose(forcing[:, 4], 1.0, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(balance[:, -3:], 0.0, rtol=0.0, atol=1e-12)
    np.testing.assert_allclose(balance[:, 13:16], -forcing[:, 1:4], rtol=0.0, atol=5e-15)

    # The wall discretization evolves the sampled parabola, exercising AB2 history.
    assert np.ptp(balance[:, 1]) > 1e-8


if __name__ == "__main__":
    test_zero_net_acceleration()
