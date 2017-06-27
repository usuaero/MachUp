"""Tests for the WingSegment geometry class"""

import json
from collections import OrderedDict
import machup.geometry as geom
import pytest
import numpy as np


@pytest.fixture
def straight_segment():
    """Returns a straight WingSegment with root at origin"""
    with open("test/geometry/testwings/wing_0.json") as file:
        wing_data = json.load(file, object_pairs_hook=OrderedDict)["wing_1"]
        wing_data["side"] = "right"
        seg = geom.WingSegment("straight_seg", wing_data)
        yield seg


WING_DIR = "test/geometry/testwings/"


# @pytest.mark.parametrize("inputfile,side,expected", [
#     (WING_DIR+"wing_0.json", "right", (0., 1., 0.)),
#     (WING_DIR+"wing_1.json", "right", (-0.173648177766692,
#                                        0.9810602621904,
#                                        -0.08583165117743)),
#     (WING_DIR+"wing_0.json", "left", (0., 1., 0.)),
#     (WING_DIR+"wing_1.json", "left", (0.173648177766692,
#                                       0.9810602621904,
#                                       0.08583165117743)),
# ])
# def test_get_unit_normal_s(inputfile, side, expected):
#     # load input file
#     with open(inputfile) as file:
#         wing_data = json.load(file, object_pairs_hook=OrderedDict)["wing_1"]
#         wing_data["side"] = side
#         seg = geom.WingSegment("seg_name", wing_data)
#
#         unit_normal = seg.get_unit_normal_s()
#         expected = np.array(expected)
#
#         assert np.allclose(unit_normal,expected,rtol=0.,atol=1e-10) is True


@pytest.mark.parametrize("inputfile,side,tip,expected", [
    (WING_DIR+"wing_1.json", "right", "right", (-0.45530792283384,
                                                3.984778792367,
                                                -1.3486229709906)),
    (WING_DIR+"wing_1.json", "right", "left", (0.25, 0.0, -1.)),
    (WING_DIR+"wing_1.json", "left", "right", (0.25, 0.0, -1.)),
    (WING_DIR+"wing_1.json", "left", "left", (-0.45530792283384,
                                              -3.984778792367,
                                              -1.3486229709906)),
])
def test_get_side_position(inputfile, side, tip, expected):
    # load input file
    with open(inputfile) as file:
        wing_data = json.load(file, object_pairs_hook=OrderedDict)["wing_1"]
        wing_data["side"] = side
        seg = geom.WingSegment("seg_name", wing_data)

        unit_normal = seg.get_side_position(tip)
        expected = np.array(expected)

        assert np.allclose(unit_normal, expected, rtol=0., atol=1e-10) is True
