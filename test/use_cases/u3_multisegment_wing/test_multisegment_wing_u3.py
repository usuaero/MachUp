"""Use Case #3: Adding a multi-segment wing to an airplane.

The following is in an example of how to add a multi-segment wing to an
Airplane object in python using the MachUp geometry module.
"""

import machup.geometry as geom
import machup.models as models
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
                                 mix={"elevator": 1.})

    main_outboard = my_plane.add_wing("main_outboard",
                                      connect_to="main_inboard",
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
    my_llmodel = models.LLModel(my_plane)

    # Solve the lifting-line model for the given condition
    controls = {
        "aileron": 4.,
        "elevator": 7.,
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
        test = np.array([-11.71264802100433,
                         -19.27362821982268,
                         -14897.23394565647,
                         -5948.612553015879,
                         -2748.571527052681,
                         62.47402578755396])

    else:
        test = np.array([-11.711041627784381,
                         -19.273700213219868,
                         -14897.070328132313,
                         -5948.6128814855092,
                         -2756.1459692244794,
                         62.370341438977121])

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-10) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-10) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-10) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-10) is True
