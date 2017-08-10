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
        self._machup_compare = False
        self._num_vortices = plane.get_num_sections()
        self._grid = LLGrid(plane, cosine_spacing)
        self._aero_data = {
            'v_loc': np.zeros((self._num_vortices, 3)),
            'u_inf': np.zeros(3),
            'rho_loc': np.ones(self._num_vortices),
            'roll_rate': 0.,
            'pitch_rate': 0.,
            'yaw_rate': 0.,
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
        self._pre_calcs = {
            "rj1i": np.zeros((self._num_vortices, self._num_vortices, 3)),
            "rj2i": np.zeros((self._num_vortices, self._num_vortices, 3)),
            "rj1i_mag": np.zeros((self._num_vortices, self._num_vortices)),
            "rj2i_mag": np.zeros((self._num_vortices, self._num_vortices)),
            "term_2": np.zeros((self._num_vortices, self._num_vortices, 3))
        }
        self._vji = None  # np.zeros((self._num_vortices,self._num_vortices,3))
        self._v_i = None  # np.zeros((self._num_vortices, 3))
        self._forces = None  # np.zeros((self._num_vortices, 3))
        self._moments = None  # np.zeros((self._num_vortices, 3))
        self._gamma = None  # np.zeros(self._num_vortices)
        self._a = None  # np.zeros((self._num_vortices, self._num_vortices))
        self._b = None  # np.zeros(self._num_vortices)
        self._pre_calculations()

    def _pre_calculations(self):
        # perform any calculations that are dependent on geometry only. This
        # allows multiple calls to the solve function to be performed faster.
        r_cp = self._grid.get_control_point_pos()
        r_1, r_2 = self._grid.get_corner_point_pos()
        u_inf = self._aero_data["u_inf"]
        rj1i = self._pre_calcs["rj1i"]
        rj2i = self._pre_calcs["rj2i"]
        rj1i_mag = self._pre_calcs["rj1i_mag"]
        rj2i_mag = self._pre_calcs["rj2i_mag"]
        term_2 = self._pre_calcs["term_2"]

        rj1i[:, :] = r_cp[:, None] - r_1
        rj2i[:, :] = r_cp[:, None] - r_2
        rj1i_mag[:, :] = np.sqrt(np.einsum('ijk,ijk->ij', rj1i, rj1i))
        rj2i_mag[:, :] = np.sqrt(np.einsum('ijk,ijk->ij', rj2i, rj2i))
        rj1irj2i_mag = np.multiply(rj1i_mag, rj2i_mag)

        n_2 = np.cross(rj1i, rj2i)*(rj1i_mag+rj2i_mag)[:, :, None]
        d_2 = rj1irj2i_mag*(rj1irj2i_mag + np.einsum('ijk,ijk->ij', rj1i, rj2i))

        with np.errstate(divide='ignore', invalid='ignore'):
            # Diagonal elements of denominator of second term should
            # all be zero but sometimes they are not due to machine
            # precision error. Setting them all to be zero guarantees
            # that the following code performs as would be expected.
            np.fill_diagonal(d_2, 0.)
            term_2[:, :] = np.true_divide(n_2, d_2[:, :, None])
            # t_2[~ np.isfinite(t_2)] = 0.
            diag = np.diag_indices(self._num_vortices)
            term_2[diag] = 0.

    def _process_aero_state(self, state):
        # Takes state data from state and constructs necessary arrays
        v_local = self._aero_data["v_loc"]
        u_inf = self._aero_data["u_inf"]
        rho_local = self._aero_data["rho_loc"]
        if state:
            if "local_state" in state:
                v_local[:] = state["local_state"][:, 1:]
                rho_local[:] = state["local_state"][:, 0]
                v_mean = np.mean(v_local, axis=0)
                u_inf[:] = v_mean/np.sqrt(v_mean[0]*v_mean[0] +
                                          v_mean[1]*v_mean[1] +
                                          v_mean[2]*v_mean[2])
            else:
                if "V_mag" not in state:
                    raise RuntimeError("Must supply 'V_mag' key and value")
                if "rho" not in state:
                    raise RuntimeError("Must supply 'rho' key and value")
                alpha = state.get("alpha", 0.)*np.pi/180.
                beta = state.get("beta", 0.)*np.pi/180.
                c_a = np.cos(alpha)
                s_a = np.sin(alpha)
                c_b = np.cos(beta)
                s_b = np.sin(beta)
                v_xyz = np.zeros(3)
                v_xyz[:] = state["V_mag"]/np.sqrt(1.-s_a*s_a*s_b*s_b)
                v_xyz[0] *= -c_a*c_b
                v_xyz[1] *= -c_a*s_b
                v_xyz[2] *= -s_a*c_b
                v_local[:] = v_xyz
                u_inf[:] = (v_xyz/np.linalg.norm(v_xyz))
                rho_local[:] = state["rho"]

            if ('roll_rate' in state or
                    'pitch_rate' in state or
                    'yaw_rate' in state):
                self._superimpose_rotation(state)

        else:
            # assume uniform flow in -x direction for now
            v_local[:, 0] = -10.
            u_inf[0] = -1.

    def _superimpose_rotation(self, state):
        # Calculate velocities due to rotation and superimpose them
        # on the freestream velocities
        r_rate = self._aero_data["roll_rate"] = state.get("roll_rate", 0.)
        p_rate = self._aero_data["pitch_rate"] = state.get("pitch_rate", 0.)
        y_rate = self._aero_data["yaw_rate"] = state.get("yaw_rate", 0.)

        rotation = np.array([r_rate, p_rate, y_rate])
        r_cp = self._grid.get_control_point_pos()
        r_cg = self._grid.get_cg_location()
        v_local = self._aero_data["v_loc"]

        r_cp_cg = r_cp - r_cg
        # v_rot = -np.cross(rotation, r_cp_cg)
        v_local -= np.cross(rotation, r_cp_cg)

    def _process_control_state(self, state):
        # Takes state data from state and constructs necessary arrays
        if state:
            delta_a = state.get("aileron", 0.)
            delta_e = state.get("elevator", 0.)
            delta_r = state.get("rudder", 0.)
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

        rj1i = self._pre_calcs["rj1i"]
        rj2i = self._pre_calcs["rj2i"]
        rj1i_mag = self._pre_calcs["rj1i_mag"]
        rj2i_mag = self._pre_calcs["rj2i_mag"]
        term_2 = self._pre_calcs["term_2"]

        term_1 = np.cross(u_inf, rj2i)
        term_1 /= (rj2i_mag*(rj2i_mag - np.inner(u_inf, rj2i)))[:, :, None]
        term_3 = np.cross(u_inf, rj1i)
        term_3 /= (rj1i_mag*(rj1i_mag - np.inner(u_inf, rj1i)))[:, :, None]

        self._vji = term_1 + term_2 - term_3
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

        v_loc_mag = np.sqrt(np.einsum('ij,ij->i', v_loc, v_loc))
        vji_un = np.einsum('ijk,ijk->ij', u_n[:, None], vji)
        self._a = (-v_loc_mag*delta_s*cl_a)[:, None]*vji_un
        v_x_dl = np.cross(v_loc, delta_l)
        np.fill_diagonal(self._a, self._a.diagonal() +
                         2.*np.sqrt(np.einsum('ij,ij->i', v_x_dl, v_x_dl)))

        # Dimensional version (no dividing by local input velocity)
        self._b = (v_loc_mag*delta_s*cl_a *
                   (np.einsum('ij,ij->i', v_loc, u_n) +
                    v_loc_mag*(delta_flap-alpha_l0)))

        # use the following b to compare with the fortran version linear solver
        if self._machup_compare:
            # pylint: disable=no-member
            # The following matches how the fortran version of machup
            # interpolates moment coefficient across each wingsegment.
            u_a = self._grid.get_unit_axial_vectors()
            cla_left = self._grid.get_left_lift_slopes()
            cla_right = self._grid.get_right_lift_slopes()
            al0_left = self._grid.get_left_zero_lift_alpha()
            al0_right = self._grid.get_right_zero_lift_alpha()
            spacing = self._grid.get_cp_spacing()
            alpha_loc = np.arctan(np.einsum('ij,ij->i', v_loc, u_n) /
                                  np.einsum('ij,ij->i', v_loc, u_a))
            cl_left = cla_left*(alpha_loc + (delta_flap-al0_left))
            cl_right = cla_right*(alpha_loc + (delta_flap-al0_right))
            cl = cl_left + spacing*(cl_right - cl_left)
            self._b = v_loc_mag*v_loc_mag*delta_s*cl

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

        self._v_i = v_loc + np.einsum('ji,ijk->ik', gamma[:, None], vji)

    def _compute_forces(self):
        # Computes the aerodynamic force at each control point.
        rho = self._aero_data["rho_loc"]
        gamma = self._gamma
        u_inf = self._aero_data["u_inf"]
        r_1, r_2 = self._grid.get_corner_point_pos()
        delta_l = r_2 - r_1
        v_i = self._v_i

        self._forces = rho[:, None]*gamma[:, None]*np.cross(v_i, delta_l)
        force_total = np.sum(self._forces, axis=0)

        drag = np.dot(force_total, u_inf)
        lift = force_total-drag*u_inf
        lift = np.sqrt(lift[0]*lift[0]+lift[1]*lift[1]+lift[2]*lift[2])

        self._results["FX"] = force_total[0]
        self._results["FY"] = force_total[1]
        self._results["FZ"] = force_total[2]
        self._results["FL"] = lift
        self._results["FD"] = drag

    def _compute_moments(self):
        # Computes the aerodynamic moment at each control point.
        rho = self._aero_data["rho_loc"]
        r_cp = self._grid.get_control_point_pos()
        u_s = self._grid.get_unit_spanwise_vectors()
        v_i = self._v_i
        cg_location = self._grid.get_cg_location()
        force = self._forces
        int_chord2 = self._grid.get_integral_chord2()
        c_m = self._local_moment_coefficient(v_i)

        v_i_mag = np.sqrt(np.einsum('ij,ij->i', v_i, v_i))
        # use the following v_i_mag to compare with fortran version
        if self._machup_compare:
            v_loc = self._aero_data["v_loc"]
            v_i_mag = np.linalg.norm(v_loc, axis=1)
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
        alpha = np.arctan(np.einsum('ij,ij->i', v_i, u_n) /
                          np.einsum('ij,ij->i', v_i, u_a))

        if self._machup_compare:
            # The following matches how the fortran version of machup
            # interpolates moment coefficient across each wingsegment.
            spacing = self._grid.get_cp_spacing()
            cma_left = self._grid.get_left_moment_slopes()
            cma_right = self._grid.get_right_moment_slopes()
            cml0_left = self._grid.get_left_zero_lift_moments()
            cml0_right = self._grid.get_right_zero_lift_moments()
            al0_left = self._grid.get_left_zero_lift_alpha()
            al0_right = self._grid.get_right_zero_lift_alpha()

            cm_left = cml0_left+cma_left*(alpha-al0_left)
            cm_right = cml0_right+cma_right*(alpha-al0_right)
            c_m = cm_left + spacing*(cm_right - cm_left) + delta_c*cm_d
        else:
            c_m = cm_l0 + cm_a*(alpha - alpha_l0) + delta_c*cm_d

        return c_m
