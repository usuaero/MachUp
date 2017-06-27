"""Tests for the Wing geometry class"""

import json
from collections import OrderedDict
import machup.geometry as geom
import pytest


def test_get_number_of_sections():
    with open("test/geometry/testwings/wing_0.json") as file:
        wing_data = json.load(file, object_pairs_hook=OrderedDict)
        for wing_key in wing_data.keys():
            wing = geom.Wing(wing_key, wing_data[wing_key])
            assert wing.get_num_sections() == 80
