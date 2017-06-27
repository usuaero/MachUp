"""Tests for the Airplane geometry class"""
# pylint: disable=redefined-outer-name

import machup.geometry as geom
import pytest


@pytest.fixture
def single_wing_plane():
    """Returns a plane for the single_wing.json example"""
    filename = "test/geometry/testairplanes/single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    return plane


def test_get_number_of_sections(single_wing_plane):
    assert single_wing_plane.get_num_sections() == 80


def test_get_wing_segments(single_wing_plane):
    wingsegments = single_wing_plane.get_wingsegments()
    assert len(wingsegments) == 2
