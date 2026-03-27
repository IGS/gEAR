from setuptools import setup, find_packages

# This is mostly here so we can install the gEAR library utilities in the same environment as the gEAR server,
# which is required for some docker setups. Ideally we can stop appending "lib" to the PYTHONPATH with this.

# The gEAR library utilities are not intended to be used outside of the gEAR server,
# so we don't need to worry about packaging them for distribution on PyPI or anything like that.

setup(
    name="gear",
    version="1.0.0",
    packages=find_packages(),
    description="gEAR library utilities",
    python_requires=">=3.10",
)