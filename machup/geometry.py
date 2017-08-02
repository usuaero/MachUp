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

Airfoil
    Defines the properties of a 2D aifoil section.

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

    def __init__(self, name="", position=None, inputfile=None):
        self.name = name
        self._wings = {}
        self._cg_loc = np.zeros(3)
        self._position = position

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
                self.name = data["plane"]["name"]
                self._cg_loc[0] = data["plane"]["CGx"]
                self._cg_loc[1] = data["plane"]["CGy"]
                self._cg_loc[2] = data["plane"]["CGz"]
                id_map = {}
                for wing_key, wing_dict in data["wings"].items():
                    connect_to = wing_dict["connect"]["ID"]
                    if connect_to == 0:
                        id_map[wing_dict["ID"]] = wing_key
                        connect = None
                        connect_loc = 'tip'
                    else:
                        connect = id_map[connect_to]
                        connect_loc = wing_dict["connect"]["location"]

                    self.addwing(wing_key,
                                 connect_to=connect,
                                 at=connect_loc,
                                 side=wing_dict["side"],
                                 delta_pos=[wing_dict["connect"]["dx"],
                                            wing_dict["connect"]["dy"],
                                            wing_dict["connect"]["dz"]],
                                 yoffset=wing_dict["connect"]["yoffset"],
                                 semispan=wing_dict["span"],
                                 sweep=wing_dict["sweep"],
                                 dihedral=wing_dict["dihedral"],
                                 mount_angle=wing_dict["mounting_angle"],
                                 washout=wing_dict["washout"],
                                 root_chord=wing_dict["root_chord"],
                                 tip_chord=wing_dict["tip_chord"],
                                 airfoils=wing_dict["airfoils"],
                                 grid=wing_dict["grid"],
                                 control=wing_dict["control"])

    def get_num_sections(self):
        """Get the total number of sections of all of the wings.

        Returns
        -------
        int
            Total number of wing sections in airplane.

        """
        total_sections = 0
        for wing in self._wings.values():
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
        for wing in self._wings.values():
            segments.extend(wing.get_wingsegments())
        return segments

    def addwing(self, name, connect_to=None, at='tip', side='both', **dims):
        """Add wing to airplane.

        Parameters
        ----------
        wing : numpy.Wing
            Wing object that describes the wing to be added.

        Returns
        -------
        None

        """
        self._wings[name] = Wing(name, side, dims)

        if connect_to:
            parent = self._wings[connect_to]
        else:
            parent = self

        self._wings[name].connect_to(parent, at)

    def get_cg_location(self):
        """Get the location of the center of gravity.

        Returns
        -------
        Numpy Array
            Cartesian coordinates of center of gravity.

        """
        return self._cg_loc

    def get_position(self, at):
        """Get the position of the Airplane in the simulation.

        This allows the user to specify a location for the plane in
        the case that it is being used in some larger simulation.

        Parameters
        ----------
        at
            Required because get_position is a common method between
            all geometry objects but it is currently ignored in the
            case of an Airplane object.

        Returns
        -------
        Position
            Position of Airplane as specified by the user.
        """
        if self._position:
            return self._position
        else:
            return np.array([0., 0., 0.])

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

    def __init__(self, wing_name, side, dims):
        self.name = wing_name
        self._side = side
        self._left_segment = None
        self._right_segment = None

        if self._side == "both":
            self._left_segment = WingSegment("left_"+self.name,
                                             "left",
                                             dims)
            self._right_segment = WingSegment("right_"+self.name,
                                              "right",
                                              dims)
        elif self._side == "right":
            self._right_segment = WingSegment("right_"+self.name,
                                              "right",
                                              dims)
        elif self._side == "left":
            self._left_segment = WingSegment("left_"+self.name,
                                             "left",
                                             dims)
        else:
            raise RuntimeError(self._side+"side specification not recognized")

    def get_num_sections(self):
        """Get the number of sections of the wings.

        Returns
        -------
        int
            Total number of wing sections in wing.

        """
        total_sections = 0
        if self._side == "both":
            total_sections += self._left_segment.get_num_sections()
            total_sections += self._right_segment.get_num_sections()
        elif self._side == "right":
            total_sections += self._right_segment.get_num_sections()
        elif self._side == "left":
            total_sections += self._left_segment.get_num_sections()

        return total_sections

    def get_wingsegments(self):
        """Get a list of all of the WingSegments in the wing.

        Returns
        -------
        list
            List of all WingSegments in wing.

        """
        if self._side == "both":
            return [self._left_segment, self._right_segment]
        elif self._side == "right":
            return [self._right_segment]
        elif self._side == "left":
            return [self._left_segment]

    def connect_to(self, parent, at):
        """Define a connection between wing and a "parent" geometry.

        This allows for a wing to be position relative to another geometry
        to facilitate grouping of geometry together. This is done for
        convenience purposes. For example, to group tail surfaces together
        so that they can all be moved by changing the position of just one
        object instead of each individual surface.

        Parameters
        ----------
        Parent
            This can be another wing or an airplane object to connect to.

        at
            This is the location to attach to on the parent object. The
            only options currently available are "tip" which signifies a
            connection to wing tips of a parent wing.
        """
        if self._side == "both":
            self._left_segment.connect_to(parent, at)
            self._right_segment.connect_to(parent, at)
        elif self._side == "left":
            self._left_segment.connect_to(parent, at)
        elif self._side == "right":
            self._right_segment.connect_to(parent, at)

    def get_position(self, at):
        """Get the position connection location (at).

        Parameters
        ----------
        at
            The point where position is being queried. The only points
            currently available are the "left_tip" and "right_tip"
            connection points.

        Returns
        -------
        Position
            Position of connection point specified by 'at'.
        """
        if at == "left_tip":
            return self._left_segment.get_position(at)
        elif at == "right_tip":
            return self._right_segment.get_position(at)
        else:
            raise RuntimeError(at+" is not a valid connection point")


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

    def __init__(self, name, side, dims):
        self.name = name
        self._parent = None
        self._connect_at = "tip"
        self._side = side
        self._delta_pos = np.array([0., 0., 0.])
        self._dimensions = {
            "yoffset": 0.,
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
        self._unpack(dims)

    def _unpack(self, dims):
        delta_pos = dims.get("delta_pos", [0., 0., 0.])
        self._delta_pos[0] = delta_pos[0]
        self._delta_pos[1] = delta_pos[1]
        self._delta_pos[2] = delta_pos[2]
        self._dimensions["yoffset"] = dims.get("yoffset", 0.)
        self._dimensions["span"] = dims["semispan"]
        self._dimensions["root_chord"] = dims["root_chord"]
        self._dimensions["tip_chord"] = dims.get("tip_chord",
                                                 dims["root_chord"])
        self._dimensions["sweep"] = dims.get("sweep", 0.)
        self._dimensions["dihedral"] = dims.get("dihedral", 0.)
        self._dimensions["mounting_angle"] = dims.get("mount_angle", 0.)
        self._dimensions["washout"] = dims.get("washout", 0.)
        self._num_sections = dims.get("grid", 40)
        airfoils = list(dims["airfoils"].keys())
        if len(airfoils) > 1:
            self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
            self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[1]])
        else:
            self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
            self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[0]])

        control_dict = dims["control"]
        if self._side == "left":
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
            if self._side == "right":
                self._control_data["mix_aileron"] = mix_dict["aileron"]
            else:
                self._control_data["mix_aileron"] = -mix_dict["aileron"]
        if "elevator" in control_dict["mix"]:
            self._control_data["mix_elevator"] = mix_dict["elevator"]
        if "rudder" in control_dict["mix"]:
            if self._side == "right":
                self._control_data["mix_rudder"] = mix_dict["rudder"]
            else:
                self._control_data["mix_rudder"] = -mix_dict["rudder"]
        self._control_data["is_sealed"] = control_dict["is_sealed"]

    def connect_to(self, parent, at):
        """Define a connection between WingSegment and a "parent" geometry.

        This allows for a WingSegment to be position relative to another
        geometry to facilitate grouping of geometry together.

        Parameters
        ----------
        Parent
            This can be another wing or an airplane object to connect to.

        at
            This is the location to attach to on the parent object. The
            only options currently available are "tip" which signifies a
            connection to wing tips of a parent wing.
        """
        self._parent = parent
        self._connect_at = self._side+"_"+at

    def get_airfoils(self):
        """Get the root and tip Airfoil of the WingSegment.

        Returns
        -------
        tuple
            The root and tip Airfoils of the segment respectively.

        """
        return self._root_airfoil, self._tip_airfoil

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

    def get_position(self, side):
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
        if side == "left_tip":
            left_pos = np.copy(self._delta_pos)

            if self._side == "left":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                left_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                left_pos[1] -= span*np.cos(dihedral)
                left_pos[1] -= self._dimensions["yoffset"]
                left_pos[2] -= span*np.sin(dihedral)
            else:
                left_pos[1] += self._dimensions["yoffset"]

            return left_pos + self._get_parent_position()
        elif side == "right_tip":
            right_pos = np.copy(self._delta_pos)

            if self._side == "right":
                span = self._dimensions["span"]
                sweep = self._dimensions["sweep"]*np.pi/180.
                dihedral = self._dimensions["dihedral"]*np.pi/180.

                right_pos[0] -= span*np.sin(sweep)/np.cos(sweep)
                right_pos[1] += span*np.cos(dihedral)
                right_pos[1] += self._dimensions["yoffset"]
                right_pos[2] -= span*np.sin(dihedral)
            else:
                right_pos[1] -= self._dimensions["yoffset"]

            return right_pos + self._get_parent_position()
        else:
            raise RuntimeError("side specified incorrectly")

    def _get_parent_position(self):
        if self._parent:
            return self._parent.get_position(self._connect_at)
        else:
            return np.array([0., 0., 0.])

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
        return self._side

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
        # self.name = airfoil_name
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
