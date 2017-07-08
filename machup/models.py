"""Implements mathmatical modules used by machup for aircraft analysis.

Each model is encapsulated in a class with a name corresponding to the
model being implemented. These classes generally take one or more
objects that define the aircraft geometry, use this data to construct
the model, and then solve using a solve() member function.

Model list
----------
LLModel
    Numerical lifting line model for rapid analysis of aircraft
    performance and stability.
"""

import numpy as np


class LLModel:
    """The numerical lifting line model.

    The LLModel implements a modern numerical lifting line algorithm
    that models the aircraft lifting surfaces as a collection of
    horseshoe vortices. This  allows for a solution of the flow field
    and the corresponding forces and moments to be rapidly obtained.
    For further explanation, please refer to references found below.

    Unlike the traditional lifting-line theory, the numerical lifting
    line algorithm allows for multiple lifting surfaces and for sweep,
    twist, and dihedral in each surface. Viscous effects can also be
    approximated through the two-dimensional airfoil data that is used
    to close the formulation.

    Parameters
    ----------
    plane : machup.plane
        Plane object which provides necessary geometry information to
        lifting line algorithm.

    Returns
    -------
    LLModel
        Returns the newly created LLModel object.

    References
    ----------
    W. F. Phillips and D. O. Snyder. "Modern Adaptation of Prandtl's
    Classic Lifting-Line Theory", Journal of Aircraft, Vol. 37, No. 4
    (2000), pp. 662-670.

    W. F. Phillips, "Flow over Multiple Lifting Surfaces," Mechanics of
    Flight, 2nd ed., Wiley, New Jersey, 2010, pp. 94 -107.

    Examples
    --------
    A simple use case is shown below

    import machup.geometry
    import machup.models

    filename = "myAirplane.json"
    myAirplane = machup.geometry.Airplane(inputfile=filename)
    model = machup.models.LLModel(myAirplane)

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

    results = model.solve(stype="linear",
                          control_state=controls,
                          aero_state=aero_state)

    """
    # pylint: disable=too-many-instance-attributes, too-few-public-methods

    def __init__(self, plane, cosine_spacing=True):
        self._num_vortices = plane.get_num_sections()
        self._grid = LLGrid(plane, cosine_spacing)
        self._aero_data = {
            'v_loc': np.zeros((self._num_vortices, 3)),
            'u_inf': np.zeros(3),
            'rho': 1.
        }
        self._control_data = {
            'delta_flap': np.zeros(self._num_vortices),
            'deflection_eff': np.zeros(self._num_vortices)
        }
        self._results = {
            "FX": 0.,
            "FY": 0.,
            "FZ": 0.,
            "FL": 0.,
            "FD": 0.,
            "l": 0.,
            "m": 0.,
            "n": 0.
        }
        self._vji = None #np.zeros((self._num_vortices, self._num_vortices, 3))
        self._v_i = None #np.zeros((self._num_vortices, 3))
        self._forces = None #np.zeros((self._num_vortices, 3))
        self._moments = None #np.zeros((self._num_vortices, 3))
        self._gamma = None #np.zeros(self._num_vortices)
        self._a = None #np.zeros((self._num_vortices, self._num_vortices))
        self._b = None #np.zeros(self._num_vortices)

    def _process_aero_state(self, state):
        # Takes state data from state and constructs necessary arrays
        v_local = self._aero_data["v_loc"]
        u_inf = self._aero_data["u_inf"]
        if state:
            v_xyz = np.zeros(3)
            c_a = np.cos(state["alpha"]*np.pi/180.)
            s_a = np.sin(state["alpha"]*np.pi/180.)
            c_b = np.cos(state["beta"]*np.pi/180.)
            s_b = np.sin(state["beta"]*np.pi/180.)
            v_xyz[:] = state["V_mag"]/np.sqrt(1.-s_a*s_a*s_b*s_b)
            v_xyz[0] *= -c_a*c_b
            v_xyz[1] *= -c_a*s_b
            v_xyz[2] *= -s_a*c_b
            v_local[:] = v_xyz
            u_inf[:] = (v_xyz/np.linalg.norm(v_xyz))

            self._aero_data["rho"] = state["rho"]

        else:
            # assume uniform flow in -x direction for now
            v_local[:, 0] = -10.
            u_inf[0] = -1.

    def _process_control_state(self, state):
        # Takes state data from state and constructs necessary arrays
        if state:
            delta_a = state["aileron"]
            delta_e = state["elevator"]
            delta_r = state["rudder"]
            mixing_a, mixing_e, mixing_r = self._grid.get_control_mix()

            self._control_data["delta_flap"] = (delta_a*mixing_a +
                                                delta_e*mixing_e +
                                                delta_r*mixing_r)*np.pi/180.

            self._compute_deflection_efficiency()

    def _compute_deflection_efficiency(self):
        # Compute flap deflection efficiency using linear fit of
        # figure 1.7.5 in Phillips
        delta_f = self._control_data["delta_flap"]
        slope = -8.71794871794872E-03
        intercept = 1.09589743589744

        df_abs = np.absolute(delta_f)*180/np.pi
        self._control_data["deflection_eff"] = np.where(df_abs > 11.,
                                                        slope*df_abs+intercept,
                                                        1.)

    # @profile
    def solve(self, stype, aero_state=None, control_state=None):
        """Solve the numerical lifting line algorithm for the provided Plane.

        Parameters
        ----------
        stype
            Specifies the type of solution. ("linear" or "nonlinear")
        aero_state : dict
            Contains angle of attack, sideslip angle, velocity
            magnitude, and air density. Dictionary keys for these are
            "alpha", "beta", "v_mag", and "rho" respectively. Note that
            the units on density must be consistent with v_mag and the
            units used in dimensioning the Plane object. For example,
            if the plane dimensions are in feet, than v_mag should be
            in ft/s and air density should be in slug/ft^3. If no
            control state is specified, than angle of attack and
            sideslip angle are assumed to be zero and v_mag is assumed
            to be 10.
        control_state : dict
            Contains aileron, elevator, and rudder deflections in
            degrees. Dictionary keys are "aileron", "elevator", and
            "rudder". If no control_state is specified than all control
            surfaces are assumed to be at zero deflection.

        Returns
        -------
        results : dict
            Python dictionary containing the resulting forces and
            moments about the X, Y, and Z axis in the standard
            body-fixed coordinate system. Dictionary keys are "FX",
            "FY", "FZ", "l", "m", and "n".

        Raises
        ------
        RuntimeError
            Raises if stype parameter is incorrectly specified or if
            solver type is not yet implemented.

        """
        self._process_aero_state(aero_state)
        self._process_control_state(control_state)

        if stype == "linear":
            self._calc_induced_velocities()
            self._setup_linear_matrices()
            self._solve_linear_system()
            self._forces_and_moments()
            return self._results
        else:
            raise RuntimeError("solver type not yet supported")

    # @profile
    def _calc_induced_velocities(self):
        # Calculates influence of each segement of each horshoe vortex
        # on each control point. See Eq. 1.9.5 in Phillip's text. Note
        # that due to this being a dimensional version, there is no
        # characteristic length used to nondimensionalize the following
        # calculations.
        # pylint: disable=too-many-locals, no-member

        r_cp = self._grid.get_control_point_pos()
        r_1, r_2 = self._grid.get_corner_point_pos()
        u_inf = self._aero_data["u_inf"]

        rj1i = r_cp[:, None] - r_1
        rj2i = r_cp[:, None] - r_2
        rj1i_mag = np.linalg.norm(rj1i, axis=2)
        rj2i_mag = np.linalg.norm(rj2i, axis=2)
        rj1irj2i_mag = np.multiply(rj1i_mag, rj2i_mag)

        n_1 = np.cross(u_inf, rj2i)
        n_2 = np.cross(rj1i, rj2i)*(rj1i_mag+rj2i_mag)[:, :, None]
        n_3 = np.cross(u_inf, rj1i)
        d_1 = rj2i_mag*(rj2i_mag - np.inner(u_inf, rj2i))
        d_2 = rj1irj2i_mag*(rj1irj2i_mag + np.sum(rj1i*rj2i, axis=2))
        d_3 = rj1i_mag*(rj1i_mag - np.inner(u_inf, rj1i))

        self._vji = n_1/d_1[:, :, None] - n_3/d_3[:, :, None]

        with np.errstate(divide='ignore', invalid='ignore'):
            # Diagonal elements of denominator of second term should
            # all be zero but sometimes they are not due to machine
            # precision error. Setting them all to be zero guarantees
            # that the following code performs as would be expected.
            np.fill_diagonal(d_2, 0.)
            t_2 = np.true_divide(n_2, d_2[:, :, None])
            t_2[~ np.isfinite(t_2)] = 0.
        self._vji += t_2

        self._vji *= 1./(4.*np.pi)

    def _setup_linear_matrices(self):
        # Builds matrices for linear solution according to a dimensional
        # version of Eq. 1.9.17 in phillips text.
        vji = self._vji
        v_loc = self._aero_data["v_loc"]

        delta_s = self._grid.get_section_areas()
        cl_a = self._grid.get_lift_slopes()
        u_n = self._grid.get_unit_normal_vectors()
        delta_flap = self._effective_flap_deflection()
        alpha_l0 = self._grid.get_zero_lift_alpha()
        r_1, r_2 = self._grid.get_corner_point_pos()
        delta_l = r_2 - r_1

        v_loc_mag = np.linalg.norm(v_loc, axis=1)
        vji_un = np.sum(u_n[:, None]*vji, axis=2)
        self._a = (-v_loc_mag*delta_s*cl_a)[:, None]*vji_un
        np.fill_diagonal(self._a, self._a.diagonal() +
                         2.*np.linalg.norm(np.cross(v_loc, delta_l), axis=1))

        # Dimensional version (no dividing by local input velocity)
        self._b = (v_loc_mag*delta_s*cl_a *
                   (np.sum(v_loc*u_n, axis=1)+v_loc_mag*(delta_flap-alpha_l0)))

        # use the following b to compare with the fortran version linear solver
        # pylint: disable=no-member
        # u_a = self._grid.get_unit_axial_vectors()
        # alpha_loc = np.arctan(np.sum(v_loc*u_n, axis=1) /
        #                       np.sum(v_loc*u_a, axis=1))
        # self._b = (v_loc_mag*v_loc_mag*delta_s*cl_a *
        #            (alpha_loc + (delta_flap-alpha_l0)))

    def _effective_flap_deflection(self):
        # computes effective flap deflection (epsilon_f*delta) as in
        # eq 1.7.12 in Phillips text.
        eps = self._grid.get_flap_effectiveness()
        def_eff = self._control_data["deflection_eff"]
        delta_c = self._control_data["delta_flap"]

        eff_deflection = eps*def_eff*delta_c

        return eff_deflection

    def _solve_linear_system(self):
        self._gamma = np.linalg.solve(self._a, self._b)

    # @profile
    def _forces_and_moments(self):
        # uses vortex strength to compute local and total forces and moments
        self._compute_velocities()
        self._compute_forces()
        self._compute_moments()

    def _compute_velocities(self):
        # Computes the total velocity at each control point using the
        # local input velocities and the local induced velocities.
        vji = self._vji
        gamma = self._gamma
        v_loc = self._aero_data["v_loc"]

        self._v_i = v_loc + np.sum(gamma[:, None]*vji, axis=1)

    def _compute_forces(self):
        # Computes the aerodynamic force at each control point.
        rho = self._aero_data["rho"]
        gamma = self._gamma
        u_inf = self._aero_data["u_inf"]
        r_1, r_2 = self._grid.get_corner_point_pos()
        delta_l = r_2 - r_1
        v_i = self._v_i

        self._forces = rho*gamma[:, None]*np.cross(v_i, delta_l)
        force_total = np.sum(self._forces, axis=0)

        drag = np.dot(force_total, u_inf)
        lift = np.linalg.norm(force_total-drag*u_inf)

        self._results["FX"] = force_total[0]
        self._results["FY"] = force_total[1]
        self._results["FZ"] = force_total[2]
        self._results["FL"] = lift
        self._results["FD"] = drag

    def _compute_moments(self):
        # Computes the aerodynamic moment at each control point.
        rho = self._aero_data["rho"]
        r_cp = self._grid.get_control_point_pos()
        u_s = self._grid.get_unit_spanwise_vectors()
        v_i = self._v_i
        cg_location = self._grid.get_cg_location()
        force = self._forces
        int_chord2 = self._integral_chord2()
        c_m = self._local_moment_coefficient(v_i)

        v_i_mag = np.linalg.norm(v_i, axis=1)
        # use the following v_i_mag to compare with fortran version
        # v_loc = self._aero_data["v_loc"]
        # v_i_mag = np.linalg.norm(v_loc, axis=1)
        dyn_pressure = -0.5*rho*v_i_mag*v_i_mag

        self._moments = (dyn_pressure*c_m*int_chord2)[:, None]*u_s
        self._moments += np.cross((r_cp - cg_location), force)
        moment_total = np.sum(self._moments, axis=0)

        self._results["l"] = moment_total[0]
        self._results["m"] = moment_total[1]
        self._results["n"] = moment_total[2]

    def _local_moment_coefficient(self, v_i):
        # computes local moment coefficient based on local velocities
        u_a = self._grid.get_unit_axial_vectors()
        u_n = self._grid.get_unit_normal_vectors()
        delta_c = self._control_data["delta_flap"]
        cm_a, cm_l0, cm_d = self._grid.get_moment_slopes()
        alpha_l0 = self._grid.get_zero_lift_alpha()

        # pylint: disable=no-member
        alpha = np.arctan(np.sum(v_i*u_n, axis=1) /
                          np.sum(v_i*u_a, axis=1))

        c_m = cm_l0 + cm_a*(alpha - alpha_l0) + delta_c*cm_d

        return c_m

    def _integral_chord2(self):
        # Computes the integral of the chord squared along the span
        r_1, r_2 = self._grid.get_corner_point_pos()
        c_1, c_2 = self._grid.get_chord_lengths()

        int_chord2 = (np.linalg.norm((r_2-r_1), axis=1) *
                      (c_2*c_2+c_1*c_2+c_1*c_1)/3.)

        return int_chord2


class LLGrid:
    """Descritizes airplane information for the lifting-line algorithm.

    Parameters
    ----------
    machup.Plane
        Plane object that contains all of the necessary information
        about aircraft geometry.

    Returns
    -------
    LLGrid
        The newly created LLGrid object.

    """

    def __init__(self, plane, cosine_spacing=True):
        # eventually num sections will be specified through LLGrid
        # and not the geometry classes
        self._plane = plane
        self._num_sections = plane.get_num_sections()
        self._wing_segments = plane.get_wingsegments()
        self._segment_slices = []
        self._uses_cosine_spacing = cosine_spacing
        self._data = {
            'r': np.zeros((self._num_sections, 3)),
            'r_1': np.zeros((self._num_sections, 3)),
            'r_2': np.zeros((self._num_sections, 3)),
            'c_1': np.zeros(self._num_sections),
            'c_2': np.zeros(self._num_sections),
            'dS': np.zeros(self._num_sections),
            'CL_a': np.zeros(self._num_sections),
            'alpha_L0': np.zeros(self._num_sections),
            'Cm_a': np.zeros(self._num_sections),
            'Cm_L0': np.zeros(self._num_sections),
            'washout': np.zeros(self._num_sections),
            'u_a': np.zeros((self._num_sections, 3)),
            'u_n': np.zeros((self._num_sections, 3)),
            'u_s': np.zeros((self._num_sections, 3)),
            'flap_eff': np.zeros(self._num_sections),
            'mixing_a': np.zeros(self._num_sections),
            'mixing_e': np.zeros(self._num_sections),
            'mixing_r': np.zeros(self._num_sections),
            'Cm_d': np.zeros(self._num_sections)}
        self._update_data()

    def _update_data(self):
        # Builds arrays using current wing list
        index = 0
        slices = self._segment_slices
        for seg in self._wing_segments:
            num_sections = seg.get_num_sections()
            cur_slice = slice(index, index+num_sections)
            slices.append(cur_slice)
            index += num_sections

            self._calc_control_points(seg, cur_slice)
            self._calc_chord(seg, cur_slice)
            self._calc_area(seg, cur_slice)
            self._calc_washout(seg, cur_slice)
            self._calc_unit_vectors(seg, cur_slice)
            self._calc_coefficients(seg, cur_slice)
            self._calc_control_surfaces(seg, cur_slice)

    def _calc_control_points(self, seg, seg_slice):
        # Builds arrays for the control point and corner positions of
        # vortices
        num_sections = seg.get_num_sections()

        if self._uses_cosine_spacing:
            cp_spacing = self._cosine_spacing(num_sections, 0.5)
            corner_spacing = self._cosine_spacing(num_sections)
        else:
            cp_spacing = self._linear_spacing(num_sections, 0.5)
            corner_spacing = self._linear_spacing(num_sections)

        self._data["r"][seg_slice] = self._calc_segment_points(seg, cp_spacing[1:])
        self._data["r_1"][seg_slice] = self._calc_segment_points(seg, corner_spacing[:-1])
        self._data["r_2"][seg_slice] = self._calc_segment_points(seg, corner_spacing[1:])

    @staticmethod
    def _cosine_spacing(num_sections, offset=0):
        # calculates the cosine spacing
        index = np.arange(num_sections+1)
        spacing = .5*(1.-np.cos((np.pi/num_sections)*(index-offset)))

        return spacing

    @staticmethod
    def _linear_spacing(num_sections, offset=0):
        # calculates the cosine spacing
        index = np.arange(num_sections+1)
        spacing = (index-offset)/num_sections

        return spacing

    def _calc_segment_points(self, seg, spacing):
        # calculates the coordinates for a spacing of points along a segment
        left_pos = seg.get_side_position("left")
        right_pos = seg.get_side_position("right")
        length = np.linalg.norm(right_pos - left_pos)
        unit_s = self._calc_unit_s(left_pos, right_pos)
        x_pos = length*spacing*unit_s[0] + left_pos[0]
        y_pos = length*spacing*unit_s[1] + left_pos[1]
        z_pos = length*spacing*unit_s[2] + left_pos[2]

        pos = np.array([x_pos, y_pos, z_pos]).T

        return pos

    @staticmethod
    def _calc_unit_s(left_tip, right_tip):
        # calculates the unit vector from the left tip along the length
        # of the wing to the right tip. Not to be confused with the unit
        # vector in the spanwise direction which doesn't take sweep into
        # account.
        u_normal = right_tip - left_tip
        u_normal = u_normal/np.linalg.norm(u_normal)

        return u_normal

    def _calc_chord(self, seg, seg_slice):
        # calculates the chord at each wing section.
        root_chord, tip_chord = seg.get_chord()
        side = seg.get_side()
        if side == "left":
            left_chord = tip_chord
            right_chord = root_chord
        else:
            left_chord = root_chord
            right_chord = tip_chord

        r_1 = self._data["r_1"][seg_slice]
        r_2 = self._data["r_2"][seg_slice]
        c_1 = self._interp_accross_segment(left_chord, right_chord, r_1[0], r_2[-1], r_1)
        c_2 = self._interp_accross_segment(left_chord, right_chord, r_1[0], r_2[-1], r_2)

        self._data["c_1"][seg_slice] = c_1
        self._data["c_2"][seg_slice] = c_2

    @staticmethod
    def _interp_accross_segment(val_left, val_right, left_pos, right_pos, points):
        # linearly interpolates values accross wingsegment and returns and array of the
        # resulting values at the given points
        distances = np.linalg.norm(points-left_pos, axis=1)
        slope = (val_right-val_left)/np.linalg.norm(right_pos-left_pos)
        values = val_left + distances*slope

        return values

    def _calc_area(self, seg, seg_slice):
        # calculates planform area of each section
        corner_1 = self._data["r_1"][seg_slice]
        corner_2 = self._data["r_2"][seg_slice]
        chord_1 = self._data["c_1"][seg_slice]
        chord_2 = self._data["c_2"][seg_slice]

        sweep = seg.get_sweep()*np.pi/180.

        spanwise_l = np.cos(sweep)*(np.linalg.norm(corner_2-corner_1, axis=1))
        self._data["dS"][seg_slice] = spanwise_l*(chord_1+chord_2)/2.

    def _calc_washout(self, seg, seg_slice):
        # calculates the linear washout along the wing
        total_washout = seg.get_washout()
        r_cp = self._data["r"][seg_slice]
        left_tip = self._data["r_1"][seg_slice][0]
        right_tip = self._data["r_2"][seg_slice][-1]

        side = seg.get_side()
        if side == "left":
            left_washout = total_washout
            right_washout = 0.
        else:
            left_washout = 0.
            right_washout = total_washout

        washout = self._interp_accross_segment(left_washout,
                                               right_washout,
                                               left_tip,
                                               right_tip,
                                               r_cp)

        self._data["washout"][seg_slice] = washout


    def _calc_unit_vectors(self, seg, seg_slice):
        # Calculates the axial, normal, and spanwise unit vectors for
        # each section.

        washout = self._data["washout"][seg_slice]
        twist = (seg.get_mounting_angle() - washout)*np.pi/180.
        dihedral = seg.get_dihedral()*np.pi/180.
        s_twist = np.sin(twist)
        c_twist = np.cos(twist)
        s_dihedral = np.sin(dihedral)
        c_dihedral = np.cos(dihedral)
        normal = np.array([-s_twist,
                           -s_dihedral*c_twist,
                           -c_dihedral*c_twist]).T
        axial = np.array([-c_twist,
                          s_dihedral*s_twist,
                          c_dihedral*s_twist]).T

        if seg.get_side() == "left":
            normal[:, 1] *= -1.
            axial[:, 1] *= -1.

        self._data["u_n"][seg_slice] = normal
        self._data["u_a"][seg_slice] = axial
        self._data["u_s"][seg_slice] = np.cross(axial, normal)

    def _calc_coefficients(self, seg, seg_slice):
        # Calculates section airfoil properties
        root_airfoil = seg.get_root_airfoil()
        lift_slope = root_airfoil.get_lift_slope()
        alpha_l0 = root_airfoil.get_zero_lift_alpha()
        moment_slope = root_airfoil.get_moment_slope()
        moment_l0 = root_airfoil.get_zero_lift_moment()

        self._data["CL_a"][seg_slice] = lift_slope
        self._data["alpha_L0"][seg_slice] = alpha_l0
        self._data["Cm_a"][seg_slice] = moment_slope
        self._data["Cm_L0"][seg_slice] = moment_l0


    def _calc_control_surfaces(self, seg, seg_slice):
        # Sets up arrays that describe control surface properties
        # I know this is pretty messy as it currently stands but it is
        # all going to get rewritten anyway in the next version.
        # pylint: disable=too-many-locals
        flap_eff = self._data['flap_eff'][seg_slice]
        mixing_a = self._data['mixing_a'][seg_slice]
        mixing_e = self._data['mixing_e'][seg_slice]
        mixing_r = self._data['mixing_r'][seg_slice]
        cm_d = self._data['Cm_d'][seg_slice]

        num_sections = seg.get_num_sections()
        if self._uses_cosine_spacing:
            spacing = self._cosine_spacing(num_sections, 0.5)[1:]
        else:
            spacing = self._linear_spacing(num_sections, 0.5)[1:]
        surf_start, surf_end = seg.get_control_surface_span()
        chord_start = seg.get_control_surface_chord()[0]
        m_a, m_e, m_r = seg.get_control_mix()
        sealed = seg.is_control_surface_sealed()

        for i in range(num_sections):
            # if control point is covered by control surface, then
            # set mixing parameters.
            if ((spacing[i] > surf_start) and
                    (spacing[i] < surf_end)):
                # pylint: disable=no-member
                mixing_a[i] = m_a
                mixing_e[i] = m_e
                mixing_r[i] = m_r
                cf_c = chord_start
                theta_f = np.arccos(2.*cf_c - 1.)
                eps_f = 1. - (theta_f - np.sin(theta_f))/np.pi
                # the following is a curve fit of Fig. 1.7.4 in
                # Phillip's book
                eta_h = 3.9598*np.arctan((cf_c+0.006527) *
                                         89.2574+4.898015)-5.18786
                if not sealed:
                    eta_h *= 0.8
                flap_eff[i] = eta_h*eps_f
                cm_d[i] = (np.sin(2.*theta_f)-2.*np.sin(theta_f))/4.

    def get_cg_location(self):
        """Get the location of the aircraft center of gravity.

        Returns
        -------
        numpy array
            X, Y, and Z coordinates of the center of gravity.

        """
        return self._plane.get_cg_location()

    def get_control_point_pos(self):
        """Get cartesian coordinates of the control points of each vortex.

        Returns
        -------
        numpy array
            Contains the cartesian coordinates of each control point.

        """
        return self._data["r"]

    def get_corner_point_pos(self):
        """Get cartesian coordinates of the corner points of each vortex.

        Returns
        -------
        tuple of numpy arrays
            Two numpy arrays that contain the Cartesian coordinates of
            the first and second corner points of each horseshoe vortex.

        """
        return self._data["r_1"], self._data["r_2"]

    def get_lift_slopes(self):
        """Get linearly interpolated lift slopes at each wing section.

        Returns
        -------
        numpy array
            Array of lift slopes at each section.

        """
        return self._data["CL_a"]

    def get_moment_slopes(self):
        """Get linearly interpolated moment slopes at each wing section.

        Returns
        -------
        Tuple
            Array of moment coefficient slopes at each section, and zero lift
            moment coefficient at each section.

        """
        return self._data["Cm_a"], self._data["Cm_L0"], self._data["Cm_d"]

    def get_chord_lengths(self):
        """Get linearly interpolated chord lengths of each wing section.

        Returns
        -------
        tuple of numpy arrays
            Two numpy arrays that contain the chord lengths at the
            first and second corners of each horseshoe vortex
            respectively.

        """
        return self._data["c_1"], self._data["c_2"]

    def get_section_areas(self):
        """Get planform areas for each wing section.

        Returns
        -------
        numpy array
            Array of planform areas at each section.

        """
        return self._data["dS"]

    def get_unit_axial_vectors(self):
        """Get unit axial vectors for each wing section.

        Returns
        -------
        numpy array
            Contains the unit axial vectors for each section.

        """
        return self._data["u_a"]

    def get_unit_normal_vectors(self):
        """Get unit normal vectors for each wing section.

        Returns
        -------
        numpy array
            Contains the unit normal vectors for each section.

        """
        return self._data["u_n"]

    def get_unit_spanwise_vectors(self):
        """Get unit spanwise vectors for each wing section.

        Returns
        -------
        numpy array
            Contains the unit spanwise vectors for each section.

        """
        return self._data["u_s"]

    def get_flap_effectiveness(self):
        """Get the section flap effectiveness for each wing section.

        Returns
        -------
        numpy array
            Array holding the section flap effectiveness for each section.

        """
        return self._data["flap_eff"]

    def get_zero_lift_alpha(self):
        """Get the zero-lift angle of attack for each wing section.

        Returns
        -------
        numpy array
            Array holding the zero-lift angle of attack for each section.

        """
        return self._data["alpha_L0"]

    def get_control_mix(self):
        """Get the control surface mixing parameters for each wing section.

        Returns
        -------
        tuple
            Tuple of arrays holding the aileron, elevator, and rudder
            mixing for each section.

        """
        m_a = self._data["mixing_a"]
        m_e = self._data["mixing_e"]
        m_r = self._data["mixing_r"]

        return m_a, m_e, m_r
