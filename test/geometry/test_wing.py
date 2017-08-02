"""Tests for the Wing geometry class"""

import json
from collections import OrderedDict
import machup.geometry as geom
import pytest


def test_get_number_of_sections():
    with open("test/geometry/testwings/wing_0.json") as file:
        wing_data = json.load(file, object_pairs_hook=OrderedDict)["wing_1"]

        dims = {
            "delta_pos": [wing_data["connect"]["dx"],
                          wing_data["connect"]["dy"],
                          wing_data["connect"]["dz"]],
            "semispan": wing_data["span"],
            "sweep": wing_data["sweep"],
            "dihedral": wing_data["dihedral"],
            "mount_angle": wing_data["mounting_angle"],
            "washout": wing_data["washout"],
            "root_chord": wing_data["root_chord"],
            "tip_chord": wing_data["tip_chord"],
            "airfoils": wing_data["airfoils"],
            "grid": wing_data["grid"],
            "control": wing_data["control"]}
        wing = geom.Wing("wing_1", "both", dims)
        assert wing.get_num_sections() == 80
