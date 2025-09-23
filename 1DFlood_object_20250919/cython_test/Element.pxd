cdef class Element: #格子のオブジェクト
    cdef double position, elev, coeff, width, depth, old_depth
    cdef object time_evo
    cdef list upnoads, dnnoads
    cdef void solve_mass_equation(self, double )
    cdef double calc_increment(self)
    cpdef double get_variable_depth(self)
    cdef double get_position(self)
    cdef double get_elev(self)
    cdef double get_width(self)
    cdef double get_coeff(self)

cdef class Runge_Kutta_4th:
    cdef int stage
    cdef double hold
    cdef list increments
    cdef void update_stage(self)
    cdef void set_hold(self,double)
    cdef void update_stage0_variables(self,Element ,double)
    cdef void update_stage1_variables(self,Element ,double)
    cdef void update_stage2_variables(self,Element ,double)
    cdef void update_stage3_variables(self,Element ,double)
    cdef void update_depth(self,Element ,double)