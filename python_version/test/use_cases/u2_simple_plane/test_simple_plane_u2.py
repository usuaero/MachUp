"""Use Case #2: Building a simple monowing Airplane.

The following is in an example of how to build a simple monowing
Airplane object in a python script or program using the MachUp
geometry module.
"""

import machup.geometry as geom
from machup import LLModel
import numpy as np


def test_u2_simple_monowing():
    myPlane = geom.Airplane()
    myPlane.cg_location(-0.29, 0., 0.25)

    # add main wing and aileron control surface
    mainWing = myPlane.add_wing("main_wing",
                                position=[0., 0., 0.],
                                semispan=4.,
                                root_chord=1.,
                                tip_chord=0.5,
                                dihedral=3.)

    mainWing.airfoil("NACA2412",
                     alpha_L0=-0.036899751,
                     CL_alpha=6.283185307,
                     Cm_L0=-0.0527,
                     Cm_alpha=-0.08,
                     CD0=0.0055,
                     CD0_L=-0.0045,
                     CD0_L2=0.01,
                     CL_max=1.4)

    mainWing.control_surface(percent_span=(0.4, 0.9),
                             percent_chord=0.25,
                             mix={"aileron": 1.})

    # add horizontal tail and elevator control surface
    horizontalTail = myPlane.add_wing("horizontal_tail",
                                      position=[-4., 0., 0.],
                                      semispan=1.,
                                      root_chord=1.,
                                      tip_chord=0.5,
                                      mount_angle=1.,
                                      sweep=20.)

    horizontalTail.airfoil("NACA0012",
                           alpha_L0=0.,
                           CL_alpha=6.283185307,
                           Cm_L0=0.,
                           Cm_alpha=0.,
                           CD0=0.0055,
                           CD0_L=-0.0045,
                           CD0_L2=0.01,
                           CL_max=1.)

    horizontalTail.control_surface(percent_span=(0., 1.),
                                   percent_chord=(0.25, 0.5),
                                   mix={"elevator": 1.})

    # add vertical tail and rudder control surface
    verticalTail = myPlane.add_wing("vertical_tail",
                                    position=[-4., 0., 0.],
                                    side="right",
                                    semispan=0.8,
                                    root_chord=1.,
                                    tip_chord=0.5,
                                    sweep=25.,
                                    dihedral=90.)

    verticalTail.airfoil("NACA0012",
                         alpha_L0=0.,
                         CL_alpha=6.283185307,
                         Cm_L0=0.,
                         Cm_alpha=0.,
                         CD0=0.0055,
                         CD0_L=-0.0045,
                         CD0_L2=0.01,
                         CL_max=1.)

    verticalTail.control_surface(percent_span=(0., 1.),
                                 percent_chord=(0.25, 0.5),
                                 mix={"rudder": 1.})

    # Generate lifting-line model for airplane
    myLLModel = LLModel(myPlane)

    # Solve the lifting-line model for the given condition
    controls = {
        "aileron": 10.,
        "elevator": 5.,
    }
    aero_state = {
        "V_mag": 100.,
        "alpha": 5.,
        "rho": 1.
    }
    results = myLLModel.solve(stype="linear",
                              control_state=controls,
                              aero_state=aero_state)

    # compare results with expected values
    if myLLModel._machup_compare:
        test = np.array([498.563823789387015,
                         -355.559967409590001,
                         -22164.650382783511304,
                         -15767.900012819352924,
                         -6069.667491323640206,
                         -485.172664364582431])
    else:
        test = np.array([488.73383601455544,
                         -355.57437639926093,
                         -22143.208875850567,
                         -15767.798773697272,
                         -6079.3470120008151,
                         -485.09007262607491])

    assert np.allclose(results["FX"], test[0], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FY"], test[1], rtol=0., atol=1e-10) is True
    assert np.allclose(results["FZ"], test[2], rtol=0., atol=1e-10) is True
    assert np.allclose(results["l"], test[3], rtol=0., atol=1e-10) is True
    assert np.allclose(results["m"], test[4], rtol=0., atol=1e-10) is True
    assert np.allclose(results["n"], test[5], rtol=0., atol=1e-10) is True
