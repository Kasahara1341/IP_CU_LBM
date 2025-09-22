import cython
import numpy as np
cdef class Node:
    cdef double q
    cdef object down_element, up_element
    cdef double get_variable_q(self)
    cdef double calc_flux(self)
    cdef void solve_momentum_equation(self)
        
