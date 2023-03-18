import os
import subprocess
import cffi

DFTB_LIB_PATH = os.environ.get('DFTB_LIB_PATH')

# external Python callback functions for cffi
PYTHON_EXTERN = """
extern "Python" void _calc_extpot_callback(void *, double *, double *);
extern "Python" void _calc_extpotgrad_callback(void *, double *, double *);
"""

# minimum requiered version of dftbplus libary
VERSION = "0.0.0"

ABORT = False

if DFTB_LIB_PATH is None:
    # install via pipy/cmake
    import pkgconfig

    library = "dftbplus"

    try:
        if not pkgconfig.exists(library):
            raise ModuleNotFoundError(f"Unable to find pkg-config package '{library}'")
        if pkgconfig.installed(library, f"< {VERSION}"):
            raise ModuleNotFoundError(
                f"Installed '{library}' version is too old, {VERSION} or"
                " newer is required!")
        kwargs = pkgconfig.parse(library)
        cflags = pkgconfig.cflags(library).split()
        # to remove the following line you have to add .h in pkg config
        header_file = os.path.join(kwargs["library_dirs"][0],
                                    "../include/dftbplus.h")
        include_header = f'#include "{header_file}"'


        # how to make sure that pkgconfig file+header are already installed?
        # 1. You have to either set "add_subdirectory(tools)" before "Coverage
        # testing related targets" in the main CMakeLists.txt file or
        # 2. use the pkg-config file from _build file as done now

    except ModuleNotFoundError:
        try:
            prefix_var = "CONDA_PREFIX"
            kwargs = dict(libraries=[library])
            cflags = []
            if prefix_var in os.environ:
                prefix = os.environ[prefix_var]
                kwargs.update(
                    include_dirs=[os.path.join(prefix, "include")],
                    library_dirs=[os.path.join(prefix, "lib")],
                    runtime_library_dirs=[os.path.join(prefix, "lib")],
                )
                cflags.append("-I" + os.path.join(prefix, "include"))
            else:
                raise ModuleNotFoundError(
                    f"Unable to find '{library}' using pkg-config or conda!")
            # one lvl down if .h in pkg-config file via cmake install
            for directory in kwargs["include_dirs"]:
                header_file = os.path.join(directory, "dftbplus.h")
                include_header = f'#include "{header_file}"'
                if os.path.exists(header_file):
                    break
            else:
                raise ModuleNotFoundError(
                    f"Unable to find headerfile of '{library}'!")
        except ModuleNotFoundError:
            ABORT = True
        # add only one include_header here

# maybe add library dir via -rpath
else:
    # build via cmake
    kwargs = dict(libraries=["dftbplus"])
    kwargs.update(
        include_dirs=[DFTB_LIB_PATH],
        library_dirs=[DFTB_LIB_PATH],
        runtime_library_dirs=[DFTB_LIB_PATH])
    cflags = []
    header_file = os.path.join(DFTB_LIB_PATH, "dftbplus.h")
    include_header = f'#include "{header_file}"'

if not ABORT:
    cc = os.environ["CC"] if "CC" in os.environ else "cc"
    p = subprocess.Popen(
        [cc, *cflags, "-E", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = p.communicate(include_header.encode())

    cdefs = out.decode()

    ffibuilder = cffi.FFI()
    ffibuilder.set_source("pythonapi.dftb_cffi", include_header, **kwargs)
    ffibuilder.cdef(cdefs + PYTHON_EXTERN)

else:
    import warnings
    warnings.warn("No header-file found! Please make sure that a header-file"
                  " can be found via pkg-config or CONDA_PREFIX. For now "
                  "the package was installed without the pythonapi! To use "
                  "the pythonapi you have to reinstall the package.")
    raise ModuleNotFoundError()
