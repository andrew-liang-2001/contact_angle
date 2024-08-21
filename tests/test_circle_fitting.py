import pytest
import numpy as np
from contact_angle.circle_fitting import fit_circle, calculate_contact_angle


@pytest.fixture
def circle_points():
    theta = np.linspace(0, np.pi, 100)  # Only generate half-circle points
    x = 4 * np.cos(theta)  # Radius is 4
    y = 3 + 4 * np.sin(theta)  # Center y is 3
    return x, y


def test_fit_circle(circle_points):
    x, y = circle_points

    xc, yc, R = fit_circle(x, y)

    # Check if the fitted parameters match the true values
    assert xc == pytest.approx(0, abs=1e-6)  # x-center should be constrained to 0
    assert yc == pytest.approx(3, abs=1e-6)
    assert R == pytest.approx(4, abs=1e-6)


def test_calculate_contact_angle():
    # Test with known values
    R = 4  # radius
    h = 4  # height from baseline to circle center

    angle = calculate_contact_angle(R, h)

    # Expected angle when h = R is 90 degrees
    assert angle == pytest.approx(90, abs=1e-6)


def test_full_process(circle_points):
    x, y = circle_points
    baseline = -1

    # First, fit the circle
    xc, yc, R = fit_circle(x, y)

    # Then calculate the contact angle
    h = yc - baseline
    angle = calculate_contact_angle(R, h)

    # The expected angle for a circle centered at (0, 3) with radius 4 and baseline at -1
    expected_angle = np.degrees(np.arccos((4 - 4) / 4))  # should be 90 degrees

    assert angle == pytest.approx(expected_angle, abs=1e-6)
