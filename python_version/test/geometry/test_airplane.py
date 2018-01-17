"""Tests for the Airplane geometry class."""
# pylint: disable=redefined-outer-name

import machup.geometry as geom
import numpy as np
import pytest


@pytest.fixture
def single_wing_plane():
    """Return a plane for the single_wing.json example."""
    filename = "test/geometry/testairplanes/single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    return plane


@pytest.fixture
def multisegment_wing_plane():
    """Return a plane for the multisegment_wing.json example."""
    filename = "test/geometry/testairplanes/multisegment_wing.json"
    plane = geom.Airplane(inputfile=filename)
    return plane


@pytest.fixture
def no_control_surface_plane():
    """Return a plane for the no_control_surace.json example."""
    filename = "test/geometry/testairplanes/no_control_surface.json"
    plane = geom.Airplane(inputfile=filename)
    return plane


def test_get_number_of_sections(single_wing_plane):
    assert single_wing_plane.get_num_sections() == 80


def test_get_wing_segments(single_wing_plane):
    wingsegments = single_wing_plane.get_wingsegments()
    assert len(wingsegments) == 2


def test_multisegment_wing(multisegment_wing_plane):
    wingsegments = multisegment_wing_plane.get_wingsegments()
    pos = {}
    for seg in wingsegments:
        pos[seg.name] = np.array([seg.get_position("left_tip"),
                                  seg.get_position("right_tip")])

    wing_root = np.array([0.5, 0., -0.1])
    first_right = np.array([-0.3499546541037,
                            3.9392310120488,
                            -0.69459271066772]) + wing_root
    first_left = np.array([-0.3499546541037,
                           -3.9392310120488,
                           -0.69459271066772]) + wing_root
    second_right = np.array([-0.3499546541037, 4., 0.]) + first_right
    second_left = np.array([-0.3499546541037, -4., 0.]) + first_left

    l_out = np.array([second_left, first_left])
    l_in = np.array([first_left, wing_root])
    r_in = np.array([wing_root, first_right])
    r_out = np.array([first_right, second_right])

    assert np.allclose(pos["left_outer"], l_out, rtol=0., atol=1e-12) is True
    assert np.allclose(pos["left_inner"], l_in, rtol=0., atol=1e-12) is True
    assert np.allclose(pos["right_inner"], r_in, rtol=0., atol=1e-12) is True
    assert np.allclose(pos["right_outer"], r_out, rtol=0., atol=1e-12) is True


def test_no_control_surf_specified(no_control_surface_plane):
    segments = no_control_surface_plane.get_wingsegments()
    for seg in segments:
        span = seg.get_control_surface_span()
        chord = seg.get_control_surface_chord()
        mix = seg.get_control_mix()
        assert span == (0., 1.)
        assert chord == (0., 0.)
