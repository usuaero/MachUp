"""Use Case #3: Adding a multi-segment wing to an airplane.

The following is in an example of how to add a multi-segment wing to an
Airplane object in python using the MachUp geometry module.
"""

import machup.geometry as geom
from machup import LLModel
import numpy as np


def test_u3_multisegment_wing():
    my_plane = geom.Airplane()

    # add multisegment wing with aileron and flap control surfaces
    main_inboard = my_plane.add_wing("main_inboard",
                                     location=[0., 0., 0.],
                                     semispan=2.,
                                     root_chord=1.,
                                     dihedral=3.)

    main_inboard.airfoil("NACA2412",
                         alpha_L0=-0.036899751,
                         CL_alpha=6.283185307,
                         Cm_L0=-0.0527,
                         Cm_alpha=-0.08,
                         CD0=0.0055,
                         CD0_L=-0.0045,
                         CD0_L2=0.01,
                         CL_max=1.4)

    main_inboard.control_surface(percent_span=[0.2, 0.9],
                                 percent_chord=0.25,
                                 mix={"flap": 1.})

    main_outboard = my_plane.add_wing("main_outboard",
                                      connect_to=("main_inboard", "tip"),
                                      semispan=2.,
                                      root_chord=1.,
                                      tip_chord=0.6,
                                      sweep=-3.)

    main_outboard.airfoil("NACA2412",
                          alpha_L0=-0.036899751,
                          CL_alpha=6.283185307,
                          Cm_L0=-0.0527,
                          Cm_alpha=-0.08,
                          CD0=0.0055,
                          CD0_L=-0.0045,
                          CD0_L2=0.01,
                          CL_max=1.4)

    main_outboard.control_surface(percent_span=[0., 0.8],
                                  percent_chord=0.25,
                                  mix={"aileron": 1.})

    # Generate lifting-line model for airplane
    my_llmodel = LLModel(my_plane)

    # Solve the lifting-line model for the given condition
    controls = {
        "aileron": 4.,
        "flap": 7.,
    }
    aero_state = {
        "V_mag": 100.,
        "alpha": 1.,
        "rho": 1.
    }
    results = my_llmodel.solve(stype="linear",
                               control_state=controls,
                               aero_state=aero_state)

    # compare results with expected values
    if my_llmodel._machup_compare:
        test = np.array([-501.548045209023599,
                         -19.359345013074183,
                         -14896.07014164796783,
                         -5942.396955478406198,
                         -2711.801354010903196,
                         55.816799463192389])
    else:
        test = np.array([-501.76214132311253,
                         -19.339235392302943,
                         -14895.830345852655,
                         -5942.3558470931075,
                         -2719.492748965411,
                         54.85153478517509])

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-10) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-10) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-10) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-10) is True
