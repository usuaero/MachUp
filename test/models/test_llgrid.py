"""Tests for the numerical lifting line grid."""
# pylint: disable=redefined-outer-name

import pytest
import numpy as np
import machup.geometry as geom
from machup import LLGrid


PLANE_DIR = "test/geometry/testairplanes/"

# Because Machup doesn't use induced velocities in Moment calculations
COMPARING_WITH_MACHUP = False


@pytest.fixture
def single_wing_grid():
    """Get a LLGrid for the single_wing.json example."""
    filename = PLANE_DIR+"single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def small_wing_grid():
    """Get a LLGrid for the straight_wing_5sect.json example."""
    filename = PLANE_DIR+"straight_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def linear_wing_grid():
    """Get LLGrid w/ linear spacing using straight_wing_5sect.json example."""
    filename = PLANE_DIR+"straight_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane, cosine_spacing=False)
    return grid


@pytest.fixture
def vertical_wing_grid():
    """Get a LLGrid for the vertical_wing_5sect.json example."""
    filename = PLANE_DIR+"vertical_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def swept_wing_grid():
    """Get a LLGrid from the swept_wing_5sect.json example."""
    filename = PLANE_DIR+"swept_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def tapered_wing_grid():
    """Get a LLGrid from the tapered_wing_5sect.json example."""
    filename = PLANE_DIR+"tapered_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def dihedral_sweep_wing_grid():
    """Get a LLGrid from the dihedral_sweep_wing.json example."""
    filename = PLANE_DIR+"dihedral_sweep_wing.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def aero_twist_wing_grid():
    """Get a LLGrid from the aero twist example."""
    filename = PLANE_DIR+"aerodynamic_twist_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


@pytest.fixture
def yoffset_wing_grid():
    """Get a LLGrid from the yoffset wing example."""
    filename = PLANE_DIR+"yoffset_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = LLGrid(plane)
    return grid


def test_integral_of_chord_squared(dihedral_sweep_wing_grid):
    controls = {
        "aileron": 10.,
        "elevator": 0.,
        "rudder": 0.
    }
    int_chord2 = dihedral_sweep_wing_grid.get_integral_chord2()

    test = np.array([3.819660112501050975e-01,
                     1.000000000000000000e+00,
                     1.236067977499789805e+00,
                     1.000000000000000000e+00,
                     3.819660112501050975e-01,
                     3.819660112501050975e-01,
                     1.000000000000000000e+00,
                     1.236067977499789805e+00,
                     1.000000000000000000e+00,
                     3.819660112501050975e-01])

    assert np.allclose(int_chord2, test, rtol=0., atol=1e-12) is True


def test_get_grid_position(single_wing_grid):
    # get vortex positions from grid
    r_pos = single_wing_grid.get_control_point_pos()
    r_1_pos, r_2_pos = single_wing_grid.get_corner_point_pos()
    # set up what control point positins should be for the single wing case
    num_sections = 40
    span = 8.
    index = np.arange(num_sections+1)
    s_cp = np.zeros((num_sections*2, 3))
    s_1 = np.zeros((num_sections*2, 3))
    s_2 = np.zeros((num_sections*2, 3))

    r_cp = (span/4.)*(1. - np.cos((np.pi/num_sections)*(index-0.5)))[1:]
    r_1 = (span/4.)*(1. - np.cos(np.pi*index/num_sections))[:-1]
    r_2 = (span/4.)*(1. - np.cos(np.pi*index/num_sections))[1:]

    s_cp[:, 1] = np.concatenate((-np.flip(r_cp, axis=0), r_cp))
    s_1[:, 1] = np.concatenate((-np.flip(r_2, axis=0), r_1))
    s_2[:, 1] = np.concatenate((-np.flip(r_1, axis=0), r_2))

    assert len(s_cp) == len(r_pos)
    assert len(s_1) == len(r_1_pos)
    assert len(s_2) == len(r_2_pos)

    assert np.allclose(r_pos, s_cp, rtol=0., atol=1e-15) is True
    assert np.allclose(r_1_pos, s_1, rtol=0., atol=1e-15) is True
    assert np.allclose(r_2_pos, s_2, rtol=0., atol=1e-15) is True


def test_get_grid_position_swept(swept_wing_grid):
    # get vortex positions from grid
    r_pos = swept_wing_grid.get_control_point_pos()
    r_1_pos, r_2_pos = swept_wing_grid.get_corner_point_pos()
    # set up what vortex positions should be for the single wing case
    num_sections = 5
    length = 5./np.cos(np.pi/4.)
    index = np.arange(num_sections+1)
    comp_cp = (length/2.)*(1. - np.cos((np.pi/num_sections)*(index-0.5)))[1:]
    comp_1 = (length/2.)*(1. - np.cos(np.pi*index/num_sections))[:-1]
    comp_2 = (length/2.)*(1. - np.cos(np.pi*index/num_sections))[1:]
    comp_cp *= np.sqrt(2.)/2.
    comp_1 *= np.sqrt(2.)/2.
    comp_2 *= np.sqrt(2.)/2.
    # comp = np.array([0.12235870926211616, 1.0305368692688170, 2.5000000000000000,
    #                  3.9694631307311816, 4.8776412907378832])
    r_cp = np.zeros((10, 3))
    r_1 = np.zeros((10, 3))
    r_2 = np.zeros((10, 3))
    r_cp[:, 0] = np.concatenate((-1.*np.flipud(comp_cp), -1.*comp_cp))
    r_cp[:, 1] = np.concatenate((-np.flipud(comp_cp), comp_cp))
    r_1[:, 0] = np.concatenate((-1.*np.flipud(comp_2), -1.*comp_1))
    r_1[:, 1] = np.concatenate((-np.flipud(comp_2), comp_1))
    r_2[:, 0] = np.concatenate((-1.*np.flipud(comp_1), -1.*comp_2))
    r_2[:, 1] = np.concatenate((-np.flipud(comp_1), comp_2))

    assert len(r_cp) == len(r_pos)
    assert len(r_1) == len(r_1_pos)
    assert len(r_2) == len(r_2_pos)

    assert np.allclose(r_pos, r_cp, rtol=0., atol=1e-15) is True
    assert np.allclose(r_1_pos, r_1, rtol=0., atol=1e-15) is True
    assert np.allclose(r_2_pos, r_2, rtol=0., atol=1e-15) is True


def test_get_grid_position_linear(linear_wing_grid):
    # get vortex positions from grid
    r_pos = linear_wing_grid.get_control_point_pos()
    r_1_pos, r_2_pos = linear_wing_grid.get_corner_point_pos()

    r_cp = np.zeros((10, 3))
    r_1 = np.zeros((10, 3))
    r_2 = np.zeros((10, 3))

    r_cp[:, 1] = np.array([-3.6, -2.8, -2., -1.2, -0.4, 0.4, 1.2, 2., 2.8, 3.6])
    r_1[:, 1] = np.array([-4., -3.2, -2.4, -1.6, -0.8, 0., 0.8, 1.6, 2.4, 3.2])
    r_2[:, 1] = np.array([-3.2, -2.4, -1.6, -0.8, 0., 0.8, 1.6, 2.4, 3.2, 4.])

    assert np.allclose(r_pos, r_cp, rtol=0., atol=1e-15) is True
    assert np.allclose(r_1_pos, r_1, rtol=0., atol=1e-15) is True
    assert np.allclose(r_2_pos, r_2, rtol=0., atol=1e-15) is True


def test_grid_linear_interp(aero_twist_wing_grid):
    lift_slopes = aero_twist_wing_grid.get_lift_slopes()

    test_vals = np.array([6.28110521889743, 6.26555720879812,
                          6.24040000000000, 6.21524279120188,
                          6.19969478110257, 6.19969478110257,
                          6.21524279120188, 6.24040000000000,
                          6.26555720879812, 6.28110521889743])

    assert np.allclose(test_vals, lift_slopes, rtol=0., atol=1e-13) is True


def test_grid_get_lift_slope(small_wing_grid):
    lift_slopes = small_wing_grid.get_lift_slopes()

    test_vals = np.zeros(10)
    test_vals[:] = 6.1976

    assert np.allclose(test_vals, lift_slopes, rtol=0., atol=1e-13) is True


def test_grid_get_area(small_wing_grid):
    area_vals = small_wing_grid.get_section_areas()

    test_vals = np.array([0.381966011250105, 1., 1.236067977499790,
                          1., 0.381966011250105, 0.381966011250105,
                          1., 1.236067977499790, 1., 0.381966011250105])

    assert np.allclose(test_vals, area_vals, rtol=0., atol=1e-13) is True


def test_grid_get_swept_area(swept_wing_grid):
    area_vals = swept_wing_grid.get_section_areas()

    test_vals = np.array([0.954915028125264, 2.500000000000000,
                          3.090169943749470, 2.500000000000000,
                          0.954915028125263, 0.954915028125263,
                          2.500000000000000, 3.090169943749470,
                          2.500000000000000, 0.954915028125264])

    assert np.allclose(test_vals, area_vals, rtol=0., atol=1e-13) is True


def test_grid_get_tapered_area(tapered_wing_grid):
    area_vals = tapered_wing_grid.get_section_areas()

    test_vals = np.array([0.200101632734447, 0.610245751406263,
                          0.927050983124842, 0.889754248593737,
                          0.372847384140710, 0.372847384140710,
                          0.889754248593737, 0.927050983124842,
                          0.610245751406263, 0.200101632734447])

    assert np.allclose(test_vals, area_vals, rtol=0., atol=1e-13) is True


def test_grid_get_tapered_chord(tapered_wing_grid):
    chord_1, chord_2 = tapered_wing_grid.get_chord_lengths()

    test_1 = np.array([0.500000000000000, 0.547745751406263,
                       0.672745751406263, 0.827254248593737,
                       0.952254248593737, 1.000000000000000,
                       0.952254248593737, 0.827254248593737,
                       0.672745751406263, 0.547745751406263])
    test_2 = np.array([0.547745751406263, 0.672745751406263,
                       0.827254248593737, 0.952254248593737,
                       1.000000000000000, 0.952254248593737,
                       0.827254248593737, 0.672745751406263,
                       0.547745751406263, 0.500000000000000])

    assert np.allclose(test_1, chord_1, rtol=0., atol=1e-13) is True
    assert np.allclose(test_2, chord_2, rtol=0., atol=1e-13) is True


def test_grid_get_normals(small_wing_grid):
    unit_a = small_wing_grid.get_unit_axial_vectors()
    unit_n = small_wing_grid.get_unit_normal_vectors()
    unit_s = small_wing_grid.get_unit_spanwise_vectors()

    test_a = np.zeros((10, 3))
    test_n = np.zeros((10, 3))
    test_s = np.zeros((10, 3))
    test_n[:, 2] = -1.
    test_a[:, 0] = -1.
    test_s[:, 1] = -1.

    assert np.allclose(test_a, unit_a, rtol=0., atol=1e-13) is True
    assert np.allclose(test_n, unit_n, rtol=0., atol=1e-13) is True
    assert np.allclose(test_s, unit_s, rtol=0., atol=1e-13) is True


def test_grid_vertical_normals(vertical_wing_grid):
    unit_a = vertical_wing_grid.get_unit_axial_vectors()
    unit_n = vertical_wing_grid.get_unit_normal_vectors()
    unit_s = vertical_wing_grid.get_unit_spanwise_vectors()

    test_a = np.zeros((5, 3))
    test_n = np.zeros((5, 3))
    test_s = np.zeros((5, 3))
    test_a[:, 0] = -1.
    test_n[:, 1] = -1.
    test_s[:, 2] = 1.

    assert np.allclose(test_a, unit_a, rtol=0., atol=1e-13) is True
    assert np.allclose(test_n, unit_n, rtol=0., atol=1e-13) is True
    assert np.allclose(test_s, unit_s, rtol=0., atol=1e-13) is True


def test_grid_yoffset(yoffset_wing_grid):
    r_cp = yoffset_wing_grid.get_control_point_pos()

    test = np.zeros((10, 3))
    test[:, 1] = [-4.402113032590307284e+00,
                  -3.675570504584945830e+00,
                  -2.500000000000000000e+00,
                  -1.324429495415053726e+00,
                  -5.978869674096929376e-01,
                  5.978869674096929376e-01,
                  1.324429495415053726e+00,
                  2.500000000000000000e+00,
                  3.675570504584945830e+00,
                  4.402113032590307284e+00]

    assert np.allclose(r_cp, test, rtol=0., atol=1e-13) is True


def test_grid_flap_effectiveness(small_wing_grid):
    flap_eff = small_wing_grid.get_flap_effectiveness()

    test_vals = np.zeros(10)
    test_vals[1] = 0.5418715959841
    test_vals[8] = 0.5418715959841

    assert np.allclose(test_vals, flap_eff, rtol=0., atol=1e-13) is True


def test_grid_get_control_mixing(small_wing_grid):
    m_a, m_e, m_r, m_f = small_wing_grid.get_control_mix()

    test_a = np.zeros(10)
    test_a[:5] = -1.
    test_a[5:] = 1.
    test_e = np.zeros(10)
    test_e[:5] = 1.
    test_e[5:] = 1.
    test_r = np.zeros(10)
    test_f = np.zeros(10)

    assert np.allclose(test_a, m_a, rtol=0., atol=1e-13) is True
    assert np.allclose(test_e, m_e, rtol=0., atol=1e-13) is True
    assert np.allclose(test_r, m_r, rtol=0., atol=1e-13) is True
    assert np.allclose(test_f, m_f, rtol=0., atol=1e-13) is True


def test_grid_get_cm_delta(small_wing_grid):
    cm_d = small_wing_grid.get_moment_slopes()[2]

    test_vals = np.zeros(10)
    test_vals[1] = -0.649519052838329
    test_vals[8] = -0.649519052838329

    assert np.allclose(test_vals, cm_d, rtol=0., atol=1e-13) is True
