from setuptools import (
    setup,
    find_packages,
    Extension
)


module = Extension(
    "fisher",
    sources=[
        "fisher/src/fisher.c",
        "fisher/src/asa159.c",
    ],
    extra_compile_args=["-shared", "-fPIC"]
)

setup(
    install_requires=[
        "numpy",
        "scipy"
    ],
    test_suite="tests",
    tests_require=["pytest"],
    packages=find_packages(),
    package_data={
      # "fisher": ["src/*.so"]
    },
    ext_modules=[module]
)
