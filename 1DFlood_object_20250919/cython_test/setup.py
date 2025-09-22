from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        # ["Element.pyx","Calc_equation.pyx","Node.pyx","RK6th.pyx"], 
        ["Element.pyx","Node.pyx","Calc_equation.pyx"], 
        language_level=3)
)
