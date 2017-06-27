"""Tests for the numerical lifting line model"""
# pylint: disable=redefined-outer-name

import pytest
import numpy as np
import machup.geometry as geom
import machup.models as mod


PLANE_DIR = "test/geometry/testairplanes/"

# Because Machup doesn't use induced velocities in Moment calculations
COMPARING_WITH_MACHUP = False


@pytest.fixture
def single_wing_model():
    """Returns an LLModel from the single_wing.json example"""
    filename = PLANE_DIR+"single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    model = mod.LLModel(plane)
    return model


@pytest.fixture
def single_wing_grid():
    """Returns a LLGrid for the single_wing.json example"""
    filename = PLANE_DIR+"single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    grid = mod.LLGrid(plane)
    return grid


@pytest.fixture
def small_wing_model():
    """Returns a LLModel from the straight_wing_5sect.json example"""
    filename = PLANE_DIR+"straight_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = mod.LLModel(plane)
    return model


@pytest.fixture
def small_wing_grid():
    """Returns a LLGrid for the straight_wing_5sect.json example"""
    filename = PLANE_DIR+"straight_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = mod.LLGrid(plane)
    return grid


@pytest.fixture
def small_plane_model():
    """Returns a LLModel from the straight_simple_plane.json example"""
    filename = PLANE_DIR+"straight_simple_plane.json"
    plane = geom.Airplane(inputfile=filename)
    model = mod.LLModel(plane)
    return model


@pytest.fixture
def vertical_wing_model():
    """Returns a LLModel from the vertical_wing_5sect.json example"""
    filename = PLANE_DIR+"vertical_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = mod.LLModel(plane)
    return model


@pytest.fixture
def vertical_wing_grid():
    """Returns a LLGrid for the vertical_wing_5sect.json example"""
    filename = PLANE_DIR+"vertical_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    grid = mod.LLGrid(plane)
    return grid


@pytest.fixture
def swept_wing_model():
    """Returns a LLModel from the swept_wing.json example"""
    filename = PLANE_DIR+"swept_wing.json"
    plane = geom.Airplane(inputfile=filename)
    model = mod.LLModel(plane)
    return model


def test_linear_solver_forces(small_wing_model):
    results = small_wing_model.solve(stype="linear")
    machup = np.array([-1.31741502082177E-003,
                       0.00000000000000E+000,
                       -1.75579916601760E-001,
                       1.75579916601760E-001,
                       1.31741502082177E-003])

    machup[:] *= 0.5*100.*8.

    assert np.allclose(results["FX"], machup[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], machup[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], machup[2], rtol=0., atol=1e-12) is True
    # assert np.allclose(results["FL"], machup[3], rtol=0., atol=1e-12) is True
    # assert np.allclose(results["FD"], machup[4], rtol=0., atol=1e-12) is True


def test_linear_solver_moments(small_wing_model):
    results = small_wing_model.solve(stype="linear")
    if COMPARING_WITH_MACHUP:
        test = np.array([0.00000000000000E+000,
                         1.86555519922892E-002,
                         8.13151629364128E-020])
        test[0] *= 0.5*100.*8.*8.
        test[1] *= 0.5*100.*8.
        test[2] *= 0.5*100.*8.*8.
    else:
        test = np.array([0.,
                         7.4601661451501435,
                         0.])

    assert np.allclose(results["l"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[2], rtol=0., atol=1e-12) is True


def test_linear_solver_elevator(small_wing_model):
    controls = {
        "aileron": 0.,
        "elevator": 5.,
        "rudder": 0.
    }
    results = small_wing_model.solve(stype="linear", control_state=controls)

    test = np.array([-3.17807038650184E-003,
                     0.00000000000000E+000,
                     -2.28270929987400E-001,
                     -3.46944695195361E-018,
                     2.54535043011844E-002,
                     1.08420217248550E-019])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[3] = 1.0658141036401503e-14
        test[4] = 1.0168895781882410e+01
        test[5] = -3.8857805861880479e-16

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_aileron(small_wing_model):
    controls = {
        "aileron": 10.,
        "elevator": 0.,
        "rudder": 0.
    }
    results = small_wing_model.solve(stype="linear", control_state=controls)

    test = np.array([-4.65134818533851E-003,
                     0.00000000000000E+000,
                     -1.75579916601760E-001,
                     -3.59111192942685E-002,
                     1.86557335259400E-002,
                     9.87942295490482E-004])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[3] = -114.9155817416592384
        test[4] = 7.4383147236614091
        test[5] = 3.1614153455695453

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_elevon(small_wing_model):
    controls = {
        "aileron": 10.,
        "elevator": 5.,
        "rudder": 0.
    }
    results = small_wing_model.solve(stype="linear", control_state=controls)

    test = np.array([-6.20212601069838E-003,
                     0.00000000000000E+000,
                     -2.25514784671844E-001,
                     -3.49719054050339E-002,
                     2.43572038409329E-002,
                     2.18754294524915E-003])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[3] = -111.9100972961084466
        test[4] = 9.6773768983579354
        test[5] = 7.0001374247972707

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_alpha(small_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 0.,
        "rho": 1.
    }
    results = small_wing_model.solve(stype="linear",
                                     aero_state=aero_state)

    test = np.array([3.66867172969953E-002,
                     0.00000000000000E+000,
                     -5.91458366575759E-001,
                     1.38777878078145E-017,
                     1.87196576972895E-001,
                     1.73472347597681E-018])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = 14.666920954528852
        test[1] = 0.0000000000000000
        test[2] = -236.371747722803264
        test[3] = 7.460698725481052e-14
        test[4] = 7.477164959091317e+01
        test[5] = -8.326672684688674e-17

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_alpha_2(small_wing_model):
    controls = {
        "aileron": 10.,
        "elevator": 5.,
        "rudder": 0.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 0.,
        "rho": 1.
    }
    results = small_wing_model.solve(stype="linear",
                                     aero_state=aero_state,
                                     control_state=controls)

    test = np.array([3.38154734344758E-002,
                     0.00000000000000E+000,
                     -6.41886571580605E-001,
                     -3.52813757440428E-002,
                     1.93096052566309E-001,
                     1.41262716500806E-003])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = 13.519613494032514
        test[1] = 0.0000000000000000
        test[2] = -256.54292560583815
        test[3] = -112.89959289203631
        test[4] = 77.027441890667745
        test[5] = 4.511154427552104

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_beta(small_wing_model):
    controls = {
        "aileron": 0.,
        "elevator": 0.,
        "rudder": 0.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 6.,
        "rho": 1.
    }
    results = small_wing_model.solve(stype="linear",
                                     aero_state=aero_state,
                                     control_state=controls)

    test = np.array([-1.33213186788693E-003,
                     0.00000000000000E+000,
                     -1.75143481969000E-001,
                     1.08411494073456E-003,
                     1.84771495138129E-002,
                     1.01266583647442E-005])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = -0.532852747154774
        test[1] = 0.0000000000000000
        test[2] = -70.057392787599895
        test[3] = 3.469167810350585
        test[4] = 7.388776042508493
        test[5] = 0.032405306767181

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_plane(small_plane_model):
    controls = {
        "aileron": 0.,
        "elevator": 0.,
        "rudder": 0.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.
    }
    results = small_plane_model.solve(stype="linear",
                                      aero_state=aero_state,
                                      control_state=controls)

    test = np.array([-7.59542201946415E-003,
                     -1.85103255176173E-019,
                     -4.17235604408375E-001,
                     -7.21842087374036E-018,
                     8.88127478064726E-002,
                     1.06248908650680E-019])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = -3.036538764064491
        test[1] = 1.190977140616451e-15
        test[2] = -1.668493398109029e+02
        test[3] = 4.272695892230263e-14
        test[4] = 3.539916934840822e+01
        test[5] = -6.736049945627343e-15

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_plane_2(small_plane_model):
    controls = {
        "aileron": -4.,
        "elevator": 2.,
        "rudder": -2.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 3.,
        "rho": 1.
    }
    results = small_plane_model.solve(stype="linear",
                                      aero_state=aero_state,
                                      control_state=controls)

    test = np.array([4.32448388743232E-002,
                     -2.33213786289756E-003,
                     -8.62494551226151E-001,
                     1.25830282639815E-002,
                     9.52224687699335E-002,
                     1.19210981043278E-003])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = 17.279692067937397
        test[1] = -0.921170869124017
        test[2] = -343.750664533646386
        test[3] = 40.226063570577125
        test[4] = 37.559760278133943
        test[5] = 3.778720609191669

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_v_stab(vertical_wing_model):
    controls = {
        "aileron": 0.,
        "elevator": 0.,
        "rudder": 0.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 5.,
        "rho": 1.
    }
    results = vertical_wing_model.solve(stype="linear",
                                        aero_state=aero_state,
                                        control_state=controls)

    test = np.array([1.14762294095641E-002,
                     -2.46361834946595E-001,
                     -1.50853116299708E-017,
                     -6.15904587366487E-002,
                     -2.29524588191281E-002,
                     -9.19861692183353E-003])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = 4.589961997081043
        test[1] = -9.84564782591163e+01
        test[2] = -6.028720547767391e-15
        test[3] = -196.912956518232647
        test[4] = -9.179923994162081
        test[5] = -29.377214061237758

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


# def test_linear_solver_sweep(swept_wing_model):
#     aero_state = {
#         "V_mag": 10.,
#         "alpha": 4.,
#         "beta": 0.,
#         "rho": 1.
#     }
#     results = swept_wing_model.solve(stype="linear",
#                                      aero_state=aero_state)
#
#     test = np.array([6.90284991426664E-003,
#                      -1.99493199737333E-017,
#                      -1.96702770434235E-001,
#                      5.55111512312578E-018,
#                      -2.53143149215631E-001,
#                      1.80411241501588E-017])
#
#     test[:] *= 0.5*100.*20.
#     test[3] *= 20.
#     test[5] *= 20.
#     if not COMPARING_WITH_MACHUP:
#         test[0] = 14.666920954528852
#         test[1] = 0.0000000000000000
#         test[2] = -236.371747722803264
#         test[3] = 7.460698725481052e-14
#         test[4] = 7.477164959091317e+01
#         test[5] = -8.326672684688674e-17
#
#     assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
#     assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
#     assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
#     assert np.allclose(results["Cl"], test[3], rtol=0., atol=1e-12) is True
#     assert np.allclose(results["Cm"], test[4], rtol=0., atol=1e-12) is True
#     assert np.allclose(results["Cn"], test[5], rtol=0., atol=1e-12) is True


def test_get_grid_position(single_wing_grid):
    # get vortex positions from grid
    r_pos = single_wing_grid.get_control_point_pos()
    r_1_pos, r_2_pos = single_wing_grid.get_corner_point_pos()
    # set up what vortex positions should be for the single wing case
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


# def test_grid_linear_interp(small_wing_grid):
#     lift_slopes = small_wing_grid.get_lift_slopes()
#
#     test_vals = np.array([6.19981469263764, 6.21625271733377, 6.24285,
#                           6.26944728266623, 6.28588530736236,
#                           6.28588530736236, 6.26944728266623, 6.24285,
#                           6.21625271733377, 6.19981469263764])
#
#     assert np.allclose(test_vals, lift_slopes, rtol=0., atol=1e-13) is True

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


def test_grid_flap_effectiveness(small_wing_grid):
    flap_eff = small_wing_grid.get_flap_effectiveness()

    test_vals = np.zeros(10)
    test_vals[1] = 0.5418715959841
    test_vals[8] = 0.5418715959841

    assert np.allclose(test_vals, flap_eff, rtol=0., atol=1e-13) is True


def test_grid_get_control_mixing(small_wing_grid):
    m_a, m_e, m_r = small_wing_grid.get_control_mix()

    test_a = np.zeros(10)
    test_a[1] = -1.
    test_a[8] = 1.
    test_e = np.zeros(10)
    test_e[1] = 1.
    test_e[8] = 1.
    test_r = np.zeros(10)

    assert np.allclose(test_a, m_a, rtol=0., atol=1e-13) is True
    assert np.allclose(test_e, m_e, rtol=0., atol=1e-13) is True
    assert np.allclose(test_r, m_r, rtol=0., atol=1e-13) is True


def test_grid_get_cm_delta(small_wing_grid):
    cm_d = small_wing_grid.get_moment_slopes()[2]

    test_vals = np.zeros(10)
    test_vals[1] = -0.649519052838329
    test_vals[8] = -0.649519052838329

    assert np.allclose(test_vals, cm_d, rtol=0., atol=1e-13) is True
