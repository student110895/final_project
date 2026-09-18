from setuptools import setup, Extension

# Define the C extension module
module = Extension(
    "symnmfmodule",
    sources=["symnmfmodule.c", "symnmf.c"],
    libraries=["m"]
)

# Execute the setup to build the module
setup(
    name="symnmfmodule",
    version="1.0",
    ext_modules=[module]
)
