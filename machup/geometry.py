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

ControlSurface
    Defines the properties of a control surface.

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

    def _buildfrominputfile(self, filename):
        # Construct airplane from information in .json file
        if not os.path.isfile(filename):
            print('Error: Connot find file "{0}". Make sure'.format(filename))
            print(' the path is correct and the file is accessible.')
            raise IOError(filename)
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

                wing = self.add_wing(wing_key,
                                     connect_to=connect,
                                     at=connect_loc,
                                     side=wing_dict["side"],
                                     position=[wing_dict["connect"]["dx"],
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
                airfoils = list(wing_dict["airfoils"])
                if len(airfoils) > 1:
                    root = wing_dict["airfoils"][airfoils[0]]["properties"]
                    tip = wing_dict["airfoils"][airfoils[1]]["properties"]
                    wing.airfoil(airfoils[0],
                                 "root",
                                 alpha_L0=root["alpha_L0"],
                                 CL_alpha=root["CL_alpha"],
                                 Cm_L0=root["Cm_L0"],
                                 Cm_alpha=root["Cm_alpha"],
                                 CD0=root["CD0"],
                                 CD0_L=root["CD0_L"],
                                 CD0_L2=root["CD0_L2"],
                                 CL_max=root["CL_max"])
                    wing.airfoil(airfoils[1],
                                 "tip",
                                 alpha_L0=tip["alpha_L0"],
                                 CL_alpha=tip["CL_alpha"],
                                 Cm_L0=tip["Cm_L0"],
                                 Cm_alpha=tip["Cm_alpha"],
                                 CD0=tip["CD0"],
                                 CD0_L=tip["CD0_L"],
                                 CD0_L2=tip["CD0_L2"],
                                 CL_max=tip["CL_max"])
                else:
                    root = wing_dict["airfoils"][airfoils[0]]["properties"]
                    wing.airfoil(airfoils[0],
                                 "both",
                                 alpha_L0=root["alpha_L0"],
                                 CL_alpha=root["CL_alpha"],
                                 Cm_L0=root["Cm_L0"],
                                 Cm_alpha=root["Cm_alpha"],
                                 CD0=root["CD0"],
                                 CD0_L=root["CD0_L"],
                                 CD0_L2=root["CD0_L2"],
                                 CL_max=root["CL_max"])

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

    def add_wing(self, name, connect_to=None, at='tip', side='both', **dims):
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

        return self._wings[name]

    def get_cg_location(self):
        """Get the location of the center of gravity.

        Returns
        -------
        Numpy Array
            Cartesian coordinates of center of gravity.

        """
        return self._cg_loc

    def cg_location(self, x_coord=None, y_coord=None, z_coord=None):
        """Set the location of the center of gravity.

        Parameters
        -------
        x_coord
            The updated x-coordinate of the center of gravity.
        y_coord
            The updated y-coordinate of the center of gravity.
        z_coord
            The updated z-coordinate of the center of gravity.

        """
        if x_coord:
            self._cg_loc[0] = x_coord
        if y_coord:
            self._cg_loc[1] = y_coord
        if z_coord:
            self._cg_loc[2] = z_coord

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
        # pylint: disable=unused-argument
        if self._position:
            return self._position

        return np.array([0., 0., 0.])


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

    def airfoil(self, name, end="both", **properties):
        """Update the properties of the root and/or tip airfoils.

        Parameters
        ----------
        end
            Can be set as 'both', 'root', or 'tip' to set both to be
            the same or modify each individually.
        name
            The name of the airfoil.

        **properties
            All of properties need to specify an Airfoil can be passed
            in as key word arguments. These include the following:
            alpha_L0 - zero-lift angle of attack,
            CL_alpha - The lift coefficient slope,
            Cm_alpha - The moment coefficient slope,
            Cm_L0 - The zero-lift coefficient of moment,
            CD0 - The zero-lift drag,
            CD0_L - The drag coefficient proportional to lift,
            CD0_L2 - The drag coefficient proportional to lift squared.
            CL_max - The max lift coefficient.

        Returns
        -------
        Airfoil
            The newly created Airfoil object.
        """
        airfoil = Airfoil(name, properties)

        if self._side == "both":
            self._left_segment.airfoil(end, airfoil)
            self._right_segment.airfoil(end, airfoil)
        elif self._side == "left":
            self._left_segment.airfoil(end, airfoil)
        elif self._side == "right":
            self._right_segment.airfoil(end, airfoil)

        return airfoil

    def control_surface(self, mix, percent_span, percent_chord, sealed=True):
        """Set the properties of the control surface.

        Parameters
        ----------
        mix : Dict
            Contains the ratio of control surface deflection to control input.
            For example, inputing {"elevator": 1.} would mix the control
            surface 100 percent with elevator control input.

        percent_span : tuple
            The spanwise location of the control surface root and tip
            respectively, normalized by semispan.

        percent_chord : tuple or float
            The control surface width at the control surface root and tip
            respectively, normalized by chord length. If a single value is
            passed in then width is assumed to be the same at both the root
            and tip.

        is_sealed : bool
            Boolean flag for whether the control surface hinge is sealed.

        Returns
        -------
        ControlSurface
            The newly created ControlSurface object.
        """

        properties = {}
        properties["mix"] = mix
        properties["span_root"] = percent_span[0]
        properties["span_tip"] = percent_span[1]
        try:
            properties["chord_root"] = percent_chord[0]
            properties["chord_tip"] = percent_chord[1]
        except TypeError:
            properties["chord_root"] = percent_chord
            properties["chord_tip"] = percent_chord
        properties["is_sealed"] = sealed

        control_surface = ControlSurface(properties)

        if self._side == "both":
            self._left_segment.control_surface(control_surface)
            self._right_segment.control_surface(control_surface)
        elif self._side == "left":
            self._left_segment.control_surface(control_surface)
        elif self._side == "right":
            self._right_segment.control_surface(control_surface)

        return control_surface


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
        # self._control_data = {
        #     "left_span": 0.,
        #     "right_span": 0.,
        #     "left_chord": 0.,
        #     "right_chord": 0.,
        #     "mix_aileron": 0.,
        #     "mix_elevator": 0.,
        #     "mix_rudder": 0.,
        #     "is_sealed": 0.
        # }
        self._root_airfoil = None
        self._tip_airfoil = None
        self._control_surface = None
        self._num_sections = 40
        self._unpack(dims)

    def _unpack(self, dims):
        delta_pos = dims.get("position", [0., 0., 0.])
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
        if "aifoils" in dims:
            airfoils = list(dims["airfoils"].keys())
            if len(airfoils) > 1:
                self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
                self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[1]])
            else:
                self._root_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
                self._tip_airfoil = Airfoil(dims["airfoils"][airfoils[0]])
        else:
            self._root_airfoil = Airfoil()
            self._tip_airfoil = Airfoil()

        if "control" in dims:
            self._control_surface = ControlSurface(dims["control"])
        else:
            self._control_surface = ControlSurface()

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

    def airfoil(self, end, airfoil):
        """Set the airfoil(s) of the wing segment.

        Parameters
        ----------
        end : str
            Whether the airfoil is being set as the root or tip airfoil.
        airfoil : Airfoil
            The airfoil to be used.
        """
        if end == "both":
            self._root_airfoil = airfoil
            self._tip_airfoil = airfoil
        elif end == "root":
            self._root_airfoil = airfoil
        elif end == "tip":
            self._tip_airfoil = airfoil

    def control_surface(self, control_surface):
        """Set the control surface for the wing segment.

        Parameters
        ----------
        ControlSurface
            The control surface to be used.
        """
        self._control_surface = control_surface

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
        if self._side == "left":
            end, start = self._control_surface.get_control_span()
            end = 1. - end
            start = 1. - start
        elif self._side == "right":
            start, end = self._control_surface.get_control_span()

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
        if self._side == "left":
            right_chord, left_chord = self._control_surface.get_control_chord()
        elif self._side == "right":
            left_chord, right_chord = self._control_surface.get_control_chord()

        return left_chord, right_chord

    def get_control_mix(self):
        """Get the mixing parameters of the WingSegment.

        Returns
        -------
        tuple
            The aileron, elevator, and rudder mixing respectively.

        """
        control_mix = self._control_surface.get_control_mix(self._side)
        mix_aileron = control_mix.get("aileron", 0.)
        mix_elevator = control_mix.get("elevator", 0.)
        mix_rudder = control_mix.get("rudder", 0.)

        return mix_aileron, mix_elevator, mix_rudder

    def is_control_surface_sealed(self):
        """Check if control surface is sealed.

        Returns
        -------
        boolean
            True if control surface has been specified as sealed.

        """
        return self._control_surface.get_is_control_sealed()


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

    def __init__(self, airfoil_name='NACA2412', airfoil_data=None):
        self.name = airfoil_name
        if airfoil_data:
            self._properties = {
                "alpha_L0": airfoil_data["alpha_L0"],
                "CL_alpha": airfoil_data["CL_alpha"],
                "CL_max": airfoil_data["CL_max"],
                "Cm_L0": airfoil_data["Cm_L0"],
                "Cm_alpha": airfoil_data["Cm_alpha"],
                "CD_0": airfoil_data["CD0"],
                "CD_L": airfoil_data["CD0_L"],
                "CD_L2": airfoil_data["CD0_L2"]
            }
        else:
            self._properties = {
                "alpha_L0": -0.0369,
                "CL_alpha": 6.2832,
                "CL_max": 1.4,
                "Cm_L0": -0.0527,
                "Cm_alpha": -0.08,
                "CD_0": 0.0055,
                "CD_L": -0.0045,
                "CD_L2": 0.01
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


class ControlSurface:
    """Defines the dimensions and properties of wing control surface.

    Parameters
    ----------

    Returns
    -------
    ControlSurface
        Returns the newly created ControlSurface object.

    Examples
    --------
    """
    def __init__(self, dimensions=None):
        self._span_root = 0.
        self._span_tip = 0.
        self._chord_root = 0.
        self._chord_tip = 0.
        self._is_sealed = True
        self._mix = {}

        if dimensions:
            self._span_root = dimensions["span_root"]
            self._span_tip = dimensions["span_tip"]
            self._chord_root = dimensions["chord_root"]
            self._chord_tip = dimensions["chord_tip"]
            self._is_sealed = dimensions["is_sealed"]

            mix_dict = dimensions["mix"]
            if "aileron" in mix_dict:
                self._mix["aileron"] = mix_dict["aileron"]
            if "elevator" in mix_dict:
                self._mix["elevator"] = mix_dict["elevator"]
            if "rudder" in mix_dict:
                self._mix["rudder"] = mix_dict["rudder"]

    def get_control_span(self):
        """Get the spanwise location of control surface,

        Returns
        -------
        tuple
            The spanwise location of the root and tip of the control
            surface, normalized by semispan.
        """
        return self._span_root, self._span_tip

    def get_control_chord(self):
        """Get the control surface width.

        Returns
        -------
        tuple
            The control surface width at the root and tip of the control
            surface, normalized by chord.
        """
        return self._chord_root, self._chord_tip

    def get_control_mix(self, side):
        """Get the control surface mixing parameters.

        Parameters
        ----------
        side
            The side of the airplane that the WingSegment is on. This
            allows the method to return the proper mixing for control
            surfaces that are deflected in an asymmetric manner.

        Returns
        -------
        control_mix : Dict
        """
        segment_mix = {}
        if "aileron" in self._mix:
            if side == "right":
                segment_mix["aileron"] = self._mix["aileron"]
            else:
                segment_mix["aileron"] = -self._mix["aileron"]
        if "elevator" in self._mix:
            segment_mix["elevator"] = self._mix["elevator"]
        if "rudder" in self._mix:
            if side == "right":
                segment_mix["rudder"] = self._mix["rudder"]
            else:
                segment_mix["rudder"] = -self._mix["rudder"]

        return segment_mix

    def get_is_control_sealed(self):
        """Check if control surface is sealed.

        Returns
        -------
        boolean
            True if control surface has been specified as sealed.

        """
        return self._is_sealed
