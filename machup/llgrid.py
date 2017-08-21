"""Implements the grid for the lifting-line model.

Sets the spacing of the horseshoe vortices accross each wing segment,
interpolates the airfoil properties accross the wing segment, and
builds the arrays that contain all of this information for the LLModel
class.

"""


import numpy as np


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
            'spacing_cp': np.zeros(self._num_sections),
            'spacing_1': np.zeros(self._num_sections),
            'spacing_2': np.zeros(self._num_sections),
            'r': np.zeros((self._num_sections, 3)),
            'r_1': np.zeros((self._num_sections, 3)),
            'r_2': np.zeros((self._num_sections, 3)),
            'c_1': np.zeros(self._num_sections),
            'c_2': np.zeros(self._num_sections),
            'dS': np.zeros(self._num_sections),
            'CL_a': np.zeros(self._num_sections),
            'CL_a_left': np.zeros(self._num_sections),
            'CL_a_right': np.zeros(self._num_sections),
            'alpha_L0': np.zeros(self._num_sections),
            'alpha_L0_left': np.zeros(self._num_sections),
            'alpha_L0_right': np.zeros(self._num_sections),
            'Cm_a': np.zeros(self._num_sections),
            'Cm_a_left': np.zeros(self._num_sections),
            'Cm_a_right': np.zeros(self._num_sections),
            'Cm_L0': np.zeros(self._num_sections),
            'Cm_L0_left': np.zeros(self._num_sections),
            'Cm_L0_right': np.zeros(self._num_sections),
            'CD0': np.zeros(self._num_sections),
            'CD0_left': np.zeros(self._num_sections),
            'CD0_right': np.zeros(self._num_sections),
            'CD1': np.zeros(self._num_sections),
            'CD1_left': np.zeros(self._num_sections),
            'CD1_right': np.zeros(self._num_sections),
            'CD2': np.zeros(self._num_sections),
            'CD2_left': np.zeros(self._num_sections),
            'CD2_right': np.zeros(self._num_sections),
            'washout': np.zeros(self._num_sections),
            'u_a': np.zeros((self._num_sections, 3)),
            'u_n': np.zeros((self._num_sections, 3)),
            'u_s': np.zeros((self._num_sections, 3)),
            'flap_eff': np.zeros(self._num_sections),
            'mixing_a': np.zeros(self._num_sections),
            'mixing_e': np.zeros(self._num_sections),
            'mixing_r': np.zeros(self._num_sections),
            'mixing_f': np.zeros(self._num_sections),
            'Cm_d': np.zeros(self._num_sections),
            'int_chord2': np.zeros(self._num_sections)}
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
        self._calc_integral_chord2()

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

        s_cp = cp_spacing[1:]
        s_1 = corner_spacing[:-1]
        s_2 = corner_spacing[1:]

        self._data["spacing_cp"][seg_slice] = s_cp
        self._data["spacing_1"][seg_slice] = s_1
        self._data["spacing_2"][seg_slice] = s_2
        self._data["r"][seg_slice] = self._calc_segment_points(seg, s_cp)
        self._data["r_1"][seg_slice] = self._calc_segment_points(seg, s_1)
        self._data["r_2"][seg_slice] = self._calc_segment_points(seg, s_2)

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
        left_pos = seg.get_position("left_tip")
        right_pos = seg.get_position("right_tip")
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

        c_1 = self._interp_accross_segment(seg_slice, left_chord, right_chord,
                                           points="corner_1")
        c_2 = self._interp_accross_segment(seg_slice, left_chord, right_chord,
                                           points="corner_2")

        self._data["c_1"][seg_slice] = c_1
        self._data["c_2"][seg_slice] = c_2

    def _interp_accross_segment(self, seg_slice, val_left, val_right,
                                points="control"):
        # linearly interpolates values accross wingsegment and returns and
        # array of the resulting values at the given points
        if points == "control":
            spacing = self._data["spacing_cp"][seg_slice]
        elif points == "corner_1":
            spacing = self._data["spacing_1"][seg_slice]
        elif points == "corner_2":
            spacing = self._data["spacing_2"][seg_slice]

        interpolated_values = val_left + spacing*(val_right - val_left)

        return interpolated_values

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

        side = seg.get_side()
        if side == "left":
            left_washout = total_washout
            right_washout = 0.
        else:
            left_washout = 0.
            right_washout = total_washout

        washout = self._interp_accross_segment(seg_slice, left_washout,
                                               right_washout)

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
        side = seg.get_side()
        if side == "right":
            left_airfoil, right_airfoil = seg.get_airfoils()
        else:
            right_airfoil, left_airfoil = seg.get_airfoils()
        cla_left = left_airfoil.get_lift_slope()
        al0_left = left_airfoil.get_zero_lift_alpha()
        cma_left = left_airfoil.get_moment_slope()
        cml0_left = left_airfoil.get_zero_lift_moment()
        cla_right = right_airfoil.get_lift_slope()
        al0_right = right_airfoil.get_zero_lift_alpha()
        cma_right = right_airfoil.get_moment_slope()
        cml0_right = right_airfoil.get_zero_lift_moment()
        cd0_left, cd1_left, cd2_left = left_airfoil.get_drag_coefficients()
        cd0_right, cd1_right, cd2_right = right_airfoil.get_drag_coefficients()

        self._data["CL_a"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                     cla_left,
                                                                     cla_right)
        self._data["CL_a_left"][seg_slice] = cla_left
        self._data["CL_a_right"][seg_slice] = cla_right

        self._data["alpha_L0"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                         al0_left,
                                                                         al0_right)
        self._data["alpha_L0_left"][seg_slice] = al0_left
        self._data["alpha_L0_right"][seg_slice] = al0_right

        self._data["Cm_a"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                     cma_left,
                                                                     cma_right)
        self._data["Cm_a_left"][seg_slice] = cma_left
        self._data["Cm_a_right"][seg_slice] = cma_right

        self._data["Cm_L0"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                      cml0_left,
                                                                      cml0_right)
        self._data["Cm_L0_left"][seg_slice] = cml0_left
        self._data["Cm_L0_right"][seg_slice] = cml0_right

        self._data["CD0"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                    cd0_left,
                                                                    cd0_right)
        self._data["CD0_left"][seg_slice] = cd0_left
        self._data["CD0_right"][seg_slice] = cd0_right

        self._data["CD1"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                    cd1_left,
                                                                    cd1_right)
        self._data["CD1_left"][seg_slice] = cd1_left
        self._data["CD1_right"][seg_slice] = cd1_right

        self._data["CD2"][seg_slice] = self._interp_accross_segment(seg_slice,
                                                                    cd2_left,
                                                                    cd2_right)
        self._data["CD2_left"][seg_slice] = cd2_left
        self._data["CD2_right"][seg_slice] = cd2_right

    def _calc_control_surfaces(self, seg, seg_slice):
        # Sets up arrays that describe control surface properties
        # I know this is pretty messy as it currently stands but it is
        # all going to get rewritten anyway in the next version.
        # pylint: disable=too-many-locals
        flap_eff = self._data['flap_eff'][seg_slice]
        cm_d = self._data['Cm_d'][seg_slice]
        spacing = self._data["spacing_cp"][seg_slice]

        num_sections = seg.get_num_sections()
        surf_start, surf_end = seg.get_control_surface_span()
        chord_start, chord_end = seg.get_control_surface_chord()
        sealed = seg.is_control_surface_sealed()

        control_chord_slope = (chord_end-chord_start)/(surf_end-surf_start)
        cf_c = (spacing-surf_start)*control_chord_slope+chord_start

        mixing_a = self._data['mixing_a'][seg_slice]
        mixing_e = self._data['mixing_e'][seg_slice]
        mixing_r = self._data['mixing_r'][seg_slice]
        mixing_f = self._data['mixing_f'][seg_slice]

        m_a, m_e, m_r, m_f = seg.get_control_mix()

        for i in range(num_sections):
            # if control point is covered by control surface, then
            # set mixing parameters.
            if ((spacing[i] > surf_start) and (spacing[i] < surf_end)):
                # pylint: disable=no-member
                # cf_c = chord_start
                theta_f = np.arccos(2.*cf_c[i] - 1.)
                eps_f = 1. - (theta_f - np.sin(theta_f))/np.pi
                # the following is a curve fit of Fig. 1.7.4 in
                # Phillip's book
                eta_h = 3.9598*np.arctan((cf_c[i]+0.006527) *
                                         89.2574+4.898015)-5.18786
                if not sealed:
                    eta_h *= 0.8
                flap_eff[i] = eta_h*eps_f
                cm_d[i] = (np.sin(2.*theta_f)-2.*np.sin(theta_f))/4.
                mixing_a[i] = m_a
                mixing_e[i] = m_e
                mixing_r[i] = m_r
                mixing_f[i] = m_f

    def _calc_integral_chord2(self):
        # Computes the integral of the chord squared along the spanwise
        # direction
        r_1, r_2 = self.get_corner_point_pos()
        u_s = self.get_unit_spanwise_vectors()
        c_1, c_2 = self.get_chord_lengths()

        delta_l = r_2-r_1
        delta_s = np.abs(np.einsum('ij,ij->i', u_s, delta_l))
        int_chord2 = (delta_s*(c_2*c_2+c_1*c_2+c_1*c_1)/3.)

        self._data["int_chord2"] = int_chord2

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

    def get_cp_spacing(self):
        """Get the spacing of the control points along the segment span.

        Returns
        -------
        numpy array

        """
        return self._data["spacing_cp"]

    def get_lift_slopes(self):
        """Get linearly interpolated lift slopes at each wing section.

        Returns
        -------
        numpy array
            Array of lift slopes at each section.

        """
        return self._data["CL_a"]

    def get_left_lift_slopes(self):
        """Get left airfoil liftslope of segment that section lays on.

        Returns
        -------
        numpy array
            Array of left airfoil lift slopes for each section.

        """
        return self._data["CL_a_left"]

    def get_right_lift_slopes(self):
        """Get right airfoil liftslope of segment that section lays on.

        Returns
        -------
        numpy array
            Array of right airfoil lift slopes for each section.

        """
        return self._data["CL_a_right"]

    def get_drag_coefficients(self):
        """Get airfoil drag coefficients of segment that section lays on.

        Returns
        -------
        numpy array
            Array of airfoil drag coefficients for each section.

        """
        return (self._data["CD0_left"],
                self._data["CD1_left"],
                self._data["CD2_left"])

    def get_right_drag_coeff(self):
        """Get right airfoil drag coefficients of segment that section lays on.

        Returns
        -------
        numpy array
            Array of right airfoil drag coefficients for each section.

        """
        return (self._data["CD0_right"],
                self._data["CD1_right"],
                self._data["CD2_right"])

    def get_left_drag_coeff(self):
        """Get left airfoil drag coefficients of segment that section lays on.

        Returns
        -------
        numpy array
            Array of left airfoil drag coefficients for each section.

        """
        return (self._data["CD0_left"],
                self._data["CD1_left"],
                self._data["CD2_left"])

    def get_moment_slopes(self):
        """Get linearly interpolated moment slopes at each wing section.

        Returns
        -------
        Tuple
            Array of moment coefficient slopes at each section, and zero lift
            moment coefficient at each section.

        """
        return self._data["Cm_a"], self._data["Cm_L0"], self._data["Cm_d"]

    def get_left_moment_slopes(self):
        """Get left airfoil moment slopes at each wing section.

        Returns
        -------
        numpy array
            Array of left airfoil moment coefficient slopes at each section.

        """
        return self._data["Cm_a_left"]

    def get_right_moment_slopes(self):
        """Get right airfoil moment slopes at each wing section.

        Returns
        -------
        numpy array
            Array of right airfoil moment coefficient slopes at each section.

        """
        return self._data["Cm_a_right"]

    def get_left_zero_lift_moments(self):
        """Get left airfoil zero-lift moments at each wing section.

        Returns
        -------
        numpy array
            Array of left airfoil zero-lift moment coefficients.

        """
        return self._data["Cm_L0_left"]

    def get_right_zero_lift_moments(self):
        """Get right airfoil zero-lift moments at each wing section.

        Returns
        -------
        numpy array
            Array of right airfoil zero-lift moment coefficients.

        """
        return self._data["Cm_L0_right"]

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

    def get_integral_chord2(self):
        """Get integral of chord squared along span for each section.

        Returns
        -------
        numpy array
            Array of resulting values for each section.

        """
        return self._data["int_chord2"]

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

    def get_left_zero_lift_alpha(self):
        """Get the left airfoil zero-lift angle of attack for each wing section.

        Returns
        -------
        numpy array
            Array holding the left airfoil zero-lift angle of attack for each
            section.

        """
        return self._data["alpha_L0_left"]

    def get_right_zero_lift_alpha(self):
        """Get the right airfoil zero-lift angle of attack for each wing section.

        Returns
        -------
        numpy array
            Array holding the right airfoil zero-lift angle of attack for each
            section.

        """
        return self._data["alpha_L0_right"]

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
        m_f = self._data["mixing_f"]

        return m_a, m_e, m_r, m_f
