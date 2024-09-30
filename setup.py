from setuptools import setup, Extension
from Cython.Build import cythonize
import subprocess
import numpy as np


# NOTE: delete me after removing GSL dependency
def get_gsl_info():
    """Use pkg-config to retrieve GSL compile and link options."""
    try:
        # Get GSL include and library paths from pkg-config
        gsl_cflags = (
            subprocess.check_output(["pkg-config", "--cflags", "gsl"])
            .decode()
            .strip()
            .split()
        )
        gsl_ldflags = (
            subprocess.check_output(["pkg-config", "--libs", "gsl"])
            .decode()
            .strip()
            .split()
        )
        return gsl_cflags, gsl_ldflags
    except subprocess.CalledProcessError:
        raise RuntimeError(
            "GSL library not found. Please install GSL on your system."
        )


gsl_cflags, gsl_ldflags = get_gsl_info()

extensions = [
    Extension(
        "solvationfe.vdos",
        sources=[
            "solvationfe/vdos_c.c",
            "solvationfe/vdos.pyx",
        ],
        include_dirs=[np.get_include()],
        extra_compile_args=[
            *gsl_cflags,
            "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION",
        ],
        extra_link_args=gsl_ldflags,
    ),
]

setup(
    name="solvationfe",
    ext_modules=cythonize(
        extensions,
        # This shouldn't be necessary but is
        include_path=["solvationfe"],
    ),
    package_dir={"solvationfe": "solvationfe"},
)
