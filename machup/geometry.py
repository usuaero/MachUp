"""Collection of classes that describe aircraft geometry.

Each geometry class groups together all of the physical information
needed to fully define a given type of geometry. This information
is used in constructing the aerodynamic models.

Geometry list
-------------
Airplane
    Defines the geometry for a single airplane.

Wing
    Defines the geometry for a wing.

WingSegment
    Defines the geometry for a straight section of wing.

"""

import os
from collections import OrderedDict
import json
import numpy as np


class Airplane:
    """Defines the geometry for an airplane.

    A plane object fully describes a single airplane and basically
    consists of a set of wings that describe any wing and tail
    surfaces and the location of the center of gravity.

    Parameters
    ----------
    name : str
        The name of the airplane object.

    inputfile : str
        The filename for a .json file input that provides all of the
        necessary information about the airplane.

    Returns
    -------
    Airplane
        Returns the newly created Airplane object.

    Raises
    ------
    RuntimeError
        If an inputfile is not specified. This is only a temporary
        measure until the primary method of construction is implemented.
    IOError
        If input filepath or filename is invalid.

    Examples
    --------

    """

    def __init__(self, name="", inputfile=None):
        self._name = name
        self._wings = []
        self._cg_loc = np.zeros(3)

        if inputfile:
            self._buildfrominputfile(inputfile)
        else:
            raise RuntimeError("Inputfile required")

    def _buildfrominputfile(self, filename):
        # Construct airplane from information in .json file
        if not os.path.isfile(filename):
            print('Error: Connot find file "{0}". Make sure'.format(filename))
            print(' the path is correct and the file is accessible.')
            raise IOError(filename)
        else:
            with open(filename) as file:
                data = json.load(file, object_pairs_hook=OrderedDict)
                self._name = data["plane"]["name"]
                self._cg_loc[0] = data["plane"]["CGx"]
                self._cg_loc[1] = data["plane"]["CGy"]
                self._cg_loc[2] = data["plane"]["CGz"]
                for wing_key in data["wings"].keys():
                    self.addwing(Wing(wing_key, data["wings"][wing_key]))

    def get_num_sections(self):
        """Get the total number of sections of all of the wings.

        Returns
        -------
        int
            Total number of wing sections in airplane.

        """
        total_sections = 0
        for wing in self._wings:
            total_sections += wing.get_num_sections()
        return total_sections

    def get_wingsegments(self):
        """Get a list of all of the WingSegments in the airplane.

        Returns
        -------
        list
            List of all WingSegments in airplane.

        """
        segments = []
        for wing in self._wings:
            segments.extend(wing.get_wingsegments())
        return segments

    def addwing(self, wing):
        """Add wing to airplane.

        Parameters
        ----------
        wing : numpy.Wing
            Wing object that describes the wing to be added.

        Returns
        -------
        None

        """
        self._wings.append(wing)

    def get_cg_location(self):
        """Get the location of the center of gravity.

        Returns
        -------
        Numpy Array
            Cartesian coordinates of center of gravity.

        """
        return self._cg_loc

    # def update_cg_loc(self, location):
    #     """Updates the location of aircraft center of gravity.

    #     Parameters
    #     ----------
    #     location : np.array

    #     Returns
    #     -------

    #     Raises
    #     ------

    #     Notes
    #     -----

    #     """


class Wing:
    """Defines the geometry for a wing.

    Note that the main wing and horizontal and vertical stabilizers can
    all be described by the same set of information and thus do not have
    seperate classes.

    Parameters
    ----------
    wing_name : str
        Name of wing to be added. This will be used if the user wants to
        access the wing object at a later time.
    wing_dict: dict
        Python dictionary that contains all of the necessary information
        to build the wing.

    Returns
    -------
    Wing
        Returns the newly created Wing object.

    Examples:
    --------

    """

    def __init__(self, wing_name, wing_dict):
        self._name = wing_name
        self._segments = []
        self._root_loc = np.array([0., 0., 0.])
        self._is_symmetric = True
        if wing_dict:
            self._buildfromdict(wing_dict)

    def _buildfromdict(self, wing_dict):
        self._root_loc[0] = wing_dict["connect"]["dx"]
        self._root_loc[1] = wing_dict["connect"]["dy"]
        self._root_loc[2] = wing_dict["connect"]["dz"]
        if wing_dict["side"] == "both":
            self._is_symmetric = True

            wing_dict["side"] = "left"
            left_wing = WingSegment("left"+self._name, wing_dict)

            wing_dict["side"] = "right"
            right_wing = WingSegment("right"+self._name, wing_dict)

            self._segments.append(left_wing)
            self._segments.append(right_wing)
        elif wing_dict["side"] == "right":
            self._is_symmetric = False

            wing_dict["side"] = "right"
            right_wing = WingSegment("right"+self._name, wing_dict)

            self._segments.append(right_wing)
        elif wing_dict["side"] == "left":
            self._is_symmetric = False

            wing_dict["side"] = "left"
            left_wing = WingSegment("left"+self._name, wing_dict)

            self._segments.append(left_wing)
        else:
            side = wing_dict["side"]
            raise RuntimeError(side+"side specification not recognized")

    def get_num_sections(self):
        """Get the number of sections of the wings.

        Returns
        -------
        int
            Total number of wing sections in wing.

        """
        total_sections = 0
        for seg in self._segments:
            total_sections += seg.get_num_sections()
        return total_sections

    def get_wingsegments(self):
        """Get a list of all of the WingSegments in the wing.

        Returns
        -------
        list
            List of all WingSegments in wing.

        """
        return self._segments


class WingSegment:
    """Defines the geometry for a wing segment.

    A wing segment is a subdivision of a wing and describes a straight
    portion of the wing. For example, a standard main wing might have two
    segments, one for the left hand portion of the wing and another for
    right hand portion. Another example might be a wing that has an
    inboard and an outboard segment each with a different sweep and
    dihedral.

    Parameters
    ----------
    name : str
        Name of WingSegment to be added. This will be used if the user wants to
        access the WingSegment object at a later time.
    wing_dict: dict
        Python dictionary that contains all of the necessary information
        to build the WingSegment.

    Returns
    -------
    WingSegment
        Returns the newly created WingSegment object.

    """

    def __init__(self, name, wing_dict):
        self._name = name
        self._root_loc = np.array([0., 0., 0.])
        self._dimensions = {
            "side": "right",
            "span": 4.,
            "root_chord": 1.,
            "tip_chord": 1.,
            "sweep": 0.,
            "dihedral": 0.,
            "mounting_angle": 0.,
            "washout": 0
        }
        self._control_data = {
            "left_span": 0.,
            "right_span": 0.,
            "left_chord": 0.,
            "right_chord": 0.,
            "mix_aileron": 0.,
            "mix_elevator": 0.,
            "mix_rudder": 0.,
            "is_sealed": 0.
        }
        self._root_airfoil = None
        self._tip_airfoil = None
        self._num_sections = 40

        self._buildfromdict(wing_dict)

    def _buildfromdict(self, wing_dict):
        self._root_loc[0] = wing_dict["connect"]["dx"]
        self._root_loc[1] = wing_dict["connect"]["dy"]
        self._root_loc[2] = wing_dict["connect"]["dz"]
        self._dimensions["side"] = wing_dict["side"]
        self._dimensions["span"] = wing_dict["span"]
        self._dimensions["root_chord"] = wing_dict["root_chord"]
        self._dimensions["tip_chord"] = wing_dict["tip_chord"]
        self._dimensions["sweep"] = wing_dict["sweep"]
        self._dimensions["dihedral"] = wing_dict["dihedral"]
        self._dimensions["mounting_angle"] = wing_dict["mounting_angle"]
        self._dimensions["washout"] = wing_dict["washout"]
        self._num_sections = wing_dict["grid"]
        airfoils = list(wing_dict["airfoils"].keys())
        if len(airfoils) > 1:
            self._root_airfoil = Airfoil(wing_dict["airfoils"][airfoils[0]])
            self._tip_airfoil = Airfoil(wing_dict["airfoils"][airfoils[1]])
        else:
            self._root_airfoil = Airfoil(wing_dict["airfoils"][airfoils[0]])
            self._tip_airfoil = Airfoil(wing_dict["airfoils"][airfoils[0]])

        control_dict = wing_dict["control"]
        if wing_dict["side"] == "left":
            self._control_data["left_span"] = 1. - control_dict["span_tip"]
            self._control_data["right_span"] = 1. - control_dict["span_root"]
            self._control_data["left_chord"] = control_dict["chord_tip"]
            self._control_data["right_chord"] = control_dict["chord_root"]
        else:
            self._control_data["left_span"] = control_dict["span_root"]
            self._control_data["right_span"] = control_dict["span_tip"]
            self._control_data["left_chord"] = control_dict["chord_root"]
            self._control_data["right_chord"] = control_dict["chord_tip"]
        mix_dict = control_dict["mix"]
        if "aileron" in control_dict["mix"]:
            if wing_dict["side"] == "right":
                self._control_data["mix_aileron"] = mix_dict["aileron"]
            else:
                self._control_data["mix_aileron"] = -mix_dict["aileron"]
        if "elevator" in control_dict["mix"]:
            self._control_data["mix_elevator"] = mix_dict["elevator"]
        if "rudder" in control_dict["mix"]:
            self._control_data["mix_rudder"] = mix_dict["rudder"]
        self._control_data["is_sealed"] = control_dict["is_sealed"]

    def get_root_airfoil(self):
        """Get the root Airfoil of the WingSegment.

        Returns
        -------
        Airfoil
            The root Airfoil of the segment.

        """
        return self._root_airfoil

    def get_num_sections(self):
        """Get the number of sections of the WingSegment.

        Returns
        -------
        int
            Total number of wing sections in segment.

        """
        return self._num_sections

    def get_span(self):
        """Get the span of the wing segment.

        Returns
        -------
        float
            The span of the wing segment.

        """
        return self._dimensions["span"]

    def get_side_position(self, side):
        """Get the position vector of a side of the WingSegment.

        Specifically, this returns the position at right or left
        wing tip at the quarter chord.

        Parameters
        ----------
        side : str
            The side of interest. ("left" or "right")

        Returns
        -------
        np.array
            The position vector of the left or right tip at the quarter chord.

        """
        if side == "left":
            left_pos = np.copy(self._root_loc)

            if self._dimensions["side"] == "left":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                left_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                left_pos[1] -= span*np.cos(dihedral)
                left_pos[2] -= span*np.sin(dihedral)

            return left_pos
        elif side == "right":
            right_pos = np.copy(self._root_loc)

            if self._dimensions["side"] == "right":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                right_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                right_pos[1] += span*np.cos(dihedral)
                right_pos[2] -= span*np.sin(dihedral)

            return right_pos
        else:
            raise RuntimeError("side specified incorrectly")

    def get_chord(self):
        """Get the root and tip chords of the wing segment.

        Returns
        -------
        tuple
            The root and tip chords of the segment respectively.

        """
        root_chord = self._dimensions["root_chord"]
        tip_chord = self._dimensions["tip_chord"]

        return root_chord, tip_chord

    def get_mounting_angle(self):
        """Get the mounting angle of the wing segment.

        Returns
        -------
        float
            The mounting angle of the segment (degrees).

        """
        return self._dimensions["mounting_angle"]

    def get_washout(self):
        """Get the washout of the wing segment.

        Returns
        -------
        float
            The washout of the segment (degrees).

        """
        return self._dimensions["washout"]

    def get_dihedral(self):
        """Get the dihedral of the wing segment.

        Returns
        -------
        float
            The dihedral of the segment (degrees).

        """
        return self._dimensions["dihedral"]

    def get_sweep(self):
        """Get the dihedral of the wing segment.

        Returns
        -------
        float
            The dihedral of the segment (degrees).

        """
        return self._dimensions["sweep"]

    def get_side(self):
        """Get the side of the airplane that the wing segment is on.

        Returns
        -------
        float
            The side of the segment.

        """
        return self._dimensions["side"]

    def get_control_surface_span(self):
        """Get location of control surface along WingSegment.

        Position of start and end of control surface are given as a
        percentage of the WingSegment.

        Returns
        -------
        tuple
            The percent span that the control surface starts at and
            ends at starting from the left side of the WingSegment.

        """
        start = self._control_data["left_span"]
        end = self._control_data["right_span"]

        return start, end

    def get_control_surface_chord(self):
        """Get the chord of the control surface.

        Chord given as a percent of the WingSegment chord.

        Returns
        -------
        tuple
            The chord on the left side of the control surface and the
            chord on the right side of the control surface.

        """
        left_control_chord = self._control_data["left_chord"]
        right_control_chord = self._control_data["right_chord"]

        return left_control_chord, right_control_chord

    def get_control_mix(self):
        """Get the mixing parameters of the WingSegment.

        Returns
        -------
        tuple
            The aileron, elevator, and rudder mixing respectively.

        """
        mix_aileron = self._control_data["mix_aileron"]
        mix_elevator = self._control_data["mix_elevator"]
        mix_rudder = self._control_data["mix_rudder"]

        return mix_aileron, mix_elevator, mix_rudder

    def is_control_surface_sealed(self):
        """Check if control surface is sealed.

        Returns
        -------
        boolean
            True if control surface has been specified as sealed.

        """
        return self._control_data["is_sealed"]


class Airfoil:
    """Defines the aerodynamic properties of an airfoil.

    Parameters
    ----------
    airfoil_data : dict
        A python dictionary that contains all fo the properties of the
        Airfoil.

    Returns
    -------
    Airfoil
        Returns the newly created Airfoil object.

    Examples:
    --------

    """

    def __init__(self, airfoil_data):  # , airfoil_name):
        # self._name = airfoil_name
        self._properties = {
            "alpha_L0": airfoil_data["properties"]["alpha_L0"],
            "CL_alpha": airfoil_data["properties"]["CL_alpha"],
            "CL_max": airfoil_data["properties"]["CL_max"],
            "Cm_L0": airfoil_data["properties"]["Cm_L0"],
            "Cm_alpha": airfoil_data["properties"]["Cm_alpha"],
            "CD_0": airfoil_data["properties"]["CD0"],
            "CD_L": airfoil_data["properties"]["CD0_L"],
            "CD_L2": airfoil_data["properties"]["CD0_L2"]
        }

    def get_lift_slope(self):
        """Get the lift coefficient slope of the Airfoil.

        Returns
        -------
        float
            The lift slope of the Airfoil (per radian).

        """
        return self._properties["CL_alpha"]

    def get_zero_lift_alpha(self):
        """Get the zero-lift angle of attack of the Airfoil.

        Returns
        -------
        float
            The zero-lift angle of attack of the Airfoil (radians).

        """
        return self._properties["alpha_L0"]

    def get_moment_slope(self):
        """Get the moment coefficient slope of the Airfoil.

        Returns
        -------
        float
            The moment slope of the Airfoil (per radians).

        """
        return self._properties["Cm_alpha"]

    def get_zero_lift_moment(self):
        """Get the zero-lift moment coefficient of the Airfoil.

        Returns
        -------
        float
            The zero-lift moment coefficient of the Airfoil.

        """
        return self._properties["Cm_L0"]
