"""Use Case #1: Read in from .json file

The following is an example of how to construct an Airplane object from
an inputfile using the Machup geometry module.

One of the simplest ways to build an Airplane object is to read in an
existing Airplane from a .json file. These files can be exported from
the web version of Machup (see aero.go.usu.edu).

Due to its handy graphical interface, one might want to construct the
Airplane in the web version and then later export the .json file to be
used in a python script that performs more complicated analysis.
"""

import machup.geometry as geom
from machup import LLModel
import numpy as np


def test_u1_read_in_from_file():
    # Generate airplane from file
    filename = "test/use_cases/u1_from_file/use_case_1.json"
    myPlane = geom.Airplane(inputfile=filename)

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
