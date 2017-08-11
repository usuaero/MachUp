"""Tests for the numerical lifting line model."""
# pylint: disable=redefined-outer-name

import pytest
import numpy as np
import machup.geometry as geom
from machup import LLModel


PLANE_DIR = "test/geometry/testairplanes/"

# Because Machup doesn't use induced velocities in Moment calculations
COMPARING_WITH_MACHUP = False


@pytest.fixture
def single_wing_model():
    """Get a LLModel from the single_wing.json example."""
    filename = PLANE_DIR+"single_wing.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def small_wing_model():
    """Get a LLModel from the straight_wing_5sect.json example."""
    filename = PLANE_DIR+"straight_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def small_plane_model():
    """Get a LLModel from the straight_simple_plane.json example."""
    filename = PLANE_DIR+"straight_simple_plane.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def vertical_wing_model():
    """Get a LLModel from the vertical_wing_5sect.json example."""
    filename = PLANE_DIR+"vertical_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def swept_wing_model():
    """Get a LLModel from the swept_wing.json example."""
    filename = PLANE_DIR+"swept_wing.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def tapered_wing_model():
    """Get a LLModel from the tapered_wing_5sect.json example."""
    filename = PLANE_DIR+"tapered_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def washout_wing_model():
    """Get a LLModel from the washout_wing_5sect.json example."""
    filename = PLANE_DIR+"washout_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def aero_twist_wing_model():
    """Get a LLModel from the aero twist example."""
    filename = PLANE_DIR+"aerodynamic_twist_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def tapered_control_surface_model():
    """Get a LLModel from the tapered controls example."""
    filename = PLANE_DIR+"tapered_control_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def yoffset_wing_model():
    """Get a LLModel from the yoffset wing example."""
    filename = PLANE_DIR+"yoffset_wing_5sect.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
    return model


@pytest.fixture
def v2_plane_model():
    """Get a LLModel from the test plane v2 example."""
    filename = PLANE_DIR+"test_plane_v2.json"
    plane = geom.Airplane(inputfile=filename)
    model = LLModel(plane)
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
                     1.73472347597681E-018,
                     5.92405147019768E-001,
                     1.50018799815696E-002])

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
        test[6] = 236.75058824974937
        test[7] = 5.9900463451106276

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FL"], test[6], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FD"], test[7], rtol=0., atol=1e-12) is True


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
                     1.01266583647442E-005,
                     1.75143537321791E-001,
                     1.324834310131E-003])

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
    assert np.allclose(results["FL"], test[6], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FD"], test[7], rtol=0., atol=1e-12) is True


def test_linear_solver_planeV1(small_plane_model):
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


def test_linear_solver_planeV1_2(small_plane_model):
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


def test_linear_solver_planeV2(v2_plane_model):
    controls = {
        "aileron": 0.,
        "elevator": 0.,
        "rudder": 0.,
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 3.,
        "rho": 1.
    }
    results = v2_plane_model.solve(stype="linear",
                                   aero_state=aero_state,
                                   control_state=controls)

    test = np.array([1.66324304012235E-002,
                     -3.02789067480966E-002,
                     -3.90971849920270E-001,
                     -6.64426829451582E-004,
                     -1.11987165095428E-002,
                     9.81477277257450E-003])

    test[:] *= 0.5*100.*430.33
    test[3] *= 37.42
    test[4] *= 11.5
    test[5] *= 37.42
    if not COMPARING_WITH_MACHUP:
        test[0] = 357.94075831126707
        test[1] = -649.63680011889164
        test[2] = -8398.5507476769344
        test[3] = -525.35227169801283
        test[4] = -2695.3836726371997
        test[5] = 7880.4532971882736

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-10) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-10) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-10) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-10) is True


def test_linear_solver_planeV2_2(v2_plane_model):
    controls = {
        "aileron": -4.,
        "elevator": 2.,
        "rudder": -2.,
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 5.,
        "beta": 3.,
        "rho": 1.25
    }
    results = v2_plane_model.solve(stype="linear",
                                   aero_state=aero_state,
                                   control_state=controls)

    test = np.array([1.59681145902187E-002,
                     -3.67421278805487E-002,
                     -3.99215952806924E-001,
                     8.29015070213223E-003,
                     -2.36727505444233E-002,
                     1.17910303897534E-002])

    if not COMPARING_WITH_MACHUP:
        test[0] = 0.015974097378399534
        test[1] = -0.036656488116542811
        test[2] = -0.39857493418178064
        test[3] = 0.0083018123459460847
        test[4] = -0.023381882886120434
        test[5] = 0.011794827810495723

    r_x = results["FX"]/(0.5*1.25*100.*430.33)
    r_y = results["FY"]/(0.5*1.25*100.*430.33)
    r_z = results["FZ"]/(0.5*1.25*100.*430.33)
    r_l = results["l"]/(0.5*1.25*100.*430.33*37.42)
    r_m = results["m"]/(0.5*1.25*100.*430.33*11.5)
    r_n = results["n"]/(0.5*1.25*100.*430.33*37.42)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


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


def test_linear_solver_sweep(swept_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 4.,
        "beta": 0.,
        "rho": 1.
    }
    results = swept_wing_model.solve(stype="linear",
                                     aero_state=aero_state)

    test = np.array([7.87330444526774E-003,
                     2.60208521396521E-018,
                     -2.11677952845681E-001,
                     -2.22044604925031E-017,
                     -2.67548444023076E-001,
                     -4.51028103753970E-018])

    # test[:] *= 0.5*100.*20.
    # test[3] *= 10.
    # test[4] *= 2.
    # test[5] *= 10.
    if not COMPARING_WITH_MACHUP:
        test[0] = 0.0078724828127284541
        test[1] = 0.0000000000000000
        test[2] = -0.21150533232311516
        test[3] = 0.0000000000000000
        test[4] = -0.2673303653732535
        test[5] = 0.0000000000000000

    r_x = results["FX"]/(0.5*100.*20.)
    r_y = results["FY"]/(0.5*100.*20.)
    r_z = results["FZ"]/(0.5*100.*20.)
    r_l = results["l"]/(0.5*100.*20.*10.)
    r_m = results["m"]/(0.5*100.*20.*2.)
    r_n = results["n"]/(0.5*100.*20.*10.)
    # rl = results["FL"]/(0.5*100.*20.)
    # rd = results["FD"]/(0.5*100.*20.)
    # print(rl, rd)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-11) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_taper(tapered_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 4.,
        "beta": 0.,
        "rho": 1.
    }
    results = tapered_wing_model.solve(stype="linear",
                                       aero_state=aero_state)

    test = np.array([2.13070222154639E-002,
                     0.00000000000000E+000,
                     -3.63577696250575E-001,
                     1.38777878078145E-017,
                     1.45431078500230E-001,
                     8.67361737988404E-019])

    if not COMPARING_WITH_MACHUP:
        test[0] = 0.021293042810169929
        test[1] = 0.0000000000000000
        test[2] = -0.36328219682066454
        test[3] = 0.0000000000000000
        test[4] = 0.14531287872826582
        test[5] = 0.0000000000000000

    r_x = results["FX"]/(0.5*100.*6.)
    r_y = results["FY"]/(0.5*100.*6.)
    r_z = results["FZ"]/(0.5*100.*6.)
    r_l = results["l"]/(0.5*100.*6.*8.)
    r_m = results["m"]/(0.5*100.*6.*1.)
    r_n = results["n"]/(0.5*100.*6.*8.)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-11) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_washout(washout_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 4.,
        "beta": 0.,
        "rho": 1.
    }
    results = washout_wing_model.solve(stype="linear",
                                       aero_state=aero_state)

    test = np.array([2.44572906233143E-002,
                     0.00000000000000E+000,
                     -5.08457565984016E-001,
                     0.00000000000000E+000,
                     1.53558495287134E-001,
                     0.00000000000000E+000])

    if not COMPARING_WITH_MACHUP:
        test[0] = 0.024450154128978562
        test[1] = 0.0000000000000000
        test[2] = -0.50818676041079702
        test[3] = 0.0000000000000000
        test[4] = 0.1534084986133617
        test[5] = 0.0000000000000000

    r_x = results["FX"]/(0.5*100.*8.)
    r_y = results["FY"]/(0.5*100.*8.)
    r_z = results["FZ"]/(0.5*100.*8.)
    r_l = results["l"]/(0.5*100.*8.*8.)
    r_m = results["m"]/(0.5*100.*8.*1.)
    r_n = results["n"]/(0.5*100.*8.*8.)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-11) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_aero_twist(aero_twist_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.
    }
    results = aero_twist_wing_model.solve(stype="linear",
                                          aero_state=aero_state)

    test = np.array([-4.36302648170593E-004,
                     0.00000000000000E+000,
                     -9.50632859560165E-002,
                     0.00000000000000E+000,
                     1.22956281796976E-002,
                     1.01643953670516E-020])

    if not COMPARING_WITH_MACHUP:
        test[0] = -0.00043889676281888524
        test[1] = 0.0000000000000000
        test[2] = -0.095479069234970571
        test[3] = 0.0000000000000000
        test[4] = 0.012264606949093567
        test[5] = 0.0000000000000000

    r_x = results["FX"]/(0.5*100.*8.)
    r_y = results["FY"]/(0.5*100.*8.)
    r_z = results["FZ"]/(0.5*100.*8.)
    r_l = results["l"]/(0.5*100.*8.*8.)
    r_m = results["m"]/(0.5*100.*8.*1.)
    r_n = results["n"]/(0.5*100.*8.*8.)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-11) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_yoffset(yoffset_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.
    }
    results = yoffset_wing_model.solve(stype="linear",
                                       aero_state=aero_state)

    test = np.array([-1.69234613071813E-003,
                     0.00000000000000E+000,
                     -1.51576434068321E-001,
                     0.00000000000000E+000,
                     8.92791526871623E-003,
                     5.42101086242752E-020])

    if not COMPARING_WITH_MACHUP:
        test[0] = -0.0016923461307181294
        test[1] = 0.0000000000000000
        test[2] = -0.15157643406832139
        test[3] = 0.0000000000000000
        test[4] = 0.0089185193285999398
        test[5] = 0.0000000000000000

    r_x = results["FX"]/(0.5*100.*8.)
    r_y = results["FY"]/(0.5*100.*8.)
    r_z = results["FZ"]/(0.5*100.*8.)
    r_l = results["l"]/(0.5*100.*8.*8.)
    r_m = results["m"]/(0.5*100.*8.*1.)
    r_n = results["n"]/(0.5*100.*8.*8.)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-11) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_linear_solver_roll_rate(small_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.,
        "roll_rate": 0.1,
        "pitch_rate": 0.0,
        "yaw_rate": 0.0
    }
    results = small_wing_model.solve(stype="linear",
                                     aero_state=aero_state)

    test = np.array([-0.00017422282851392923,
                     0.00000000000000E+000,
                     -0.17562015983742105,
                     -0.023501253162496294,
                     0.018643763169761726,
                     -0.00035428728093468378])

    if not COMPARING_WITH_MACHUP:
        test[0] = -0.00017425831108772841
        test[1] = 0.00000000000000E+000
        test[2] = -0.17562015983742105
        test[3] = -0.023497862298498234
        test[4] = 0.018656810096338853
        test[5] = -0.00035439363246730227

    r_x = results["FX"]/(0.5*100.*8.)
    r_y = results["FY"]/(0.5*100.*8.)
    r_z = results["FZ"]/(0.5*100.*8.)
    r_l = results["l"]/(0.5*100.*8.*8.)
    r_m = results["m"]/(0.5*100.*8.*1.)
    r_n = results["n"]/(0.5*100.*8.*8.)

    assert np.allclose(r_x, test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(r_y, test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(r_z, test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(r_l, test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(r_m, test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(r_n, test[5], rtol=0., atol=1e-12) is True


def test_roll_rate_velocities(small_wing_model):
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.,
        "roll_rate": 0.1,
        "pitch_rate": 0.0,
        "yaw_rate": 0.0
    }
    results = small_wing_model.solve(stype="linear",
                                     aero_state=aero_state)
    v_local = small_wing_model._aero_data["v_loc"]

    test = np.zeros((10, 3))
    test[:, 0] = -10.
    test[:, 2] = np.array([0.390211303259031,
                           0.317557050458495,
                           0.200000000000000,
                           0.082442949541505,
                           0.009788696740969,
                           -0.009788696740969,
                           -0.082442949541505,
                           -0.200000000000000,
                           -0.317557050458495,
                           -0.390211303259031])

    assert np.allclose(v_local, test, rtol=0., atol=1e-14) is True


def test_linear_solver_tapered_controls(tapered_control_surface_model):
    controls = {
        "aileron": 5.,
        "elevator": 0.,
        "rudder": 0.
    }
    aero_state = {
        "V_mag": 10.,
        "alpha": 0.,
        "beta": 0.,
        "rho": 1.
    }
    results = tapered_control_surface_model.solve(stype="linear",
                                                  aero_state=aero_state,
                                                  control_state=controls)

    test = np.array([-3.93879262490376E-003,
                     0.00000000000000E+000,
                     -1.75579916601760E-001,
                     -3.93486689403467E-002,
                     1.86556083110728E-002,
                     8.93138706195071E-004])

    test[:] *= 0.5*100.*8.
    test[3] *= 8.
    test[5] *= 8.
    if not COMPARING_WITH_MACHUP:
        test[0] = -1.5755170499615059
        test[1] = 0.
        test[2] = -70.231966640704186
        test[3] = -125.91574060910946
        test[4] = 7.450550489412648
        test[5] = 2.8580438598242295

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-12) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-12) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-12) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-12) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-12) is True


def test_local_velocity_input(small_wing_model):
    local_state = np.array([[1.00, -10.2, 0.00, -0.30],
                           [1.05, -10.6, 0.02, -0.32],
                           [1.12, -10.1, 0.02, -0.32],
                           [1.30, -10.2, 0.04, -0.33],
                           [1.01, -10.3, 0.02, -0.31],
                           [1.04, -10.3, 0.02, -0.31],
                           [1.02, -10.6, 0.01, -0.31],
                           [1.05, -10.6, 0.05, -0.34],
                           [1.15, -10.3, 0.02, -0.32],
                           [1.05, -10.4, 0.02, -0.32]])

    aero_state = {
        "local_state": local_state
    }
    small_wing_model.solve(stype='linear', aero_state=aero_state)

    v_loc = small_wing_model._aero_data["v_loc"]
    rho_loc = small_wing_model._aero_data["rho_loc"]
    u_mean = small_wing_model._aero_data["u_inf"]
    u_test = [-0.999526990162502, 0.002122547662507, -0.03068046166715]

    assert np.allclose(v_loc, local_state[:, 1:], rtol=0., atol=1e-12) is True
    assert np.allclose(rho_loc, local_state[:, 0], rtol=0., atol=1e-12) is True
    assert np.allclose(u_mean, u_test, rtol=0., atol=1e-12) is True
