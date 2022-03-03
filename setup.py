from setuptools import setup, find_packages


setup(
    install_requires=[
        "numpy",
        "scipy"
    ],
    test_suite="tests",
    tests_require=["pytest"],
    packages=find_packages(),
    package_data={
      "fisher": ["src/*.so"]
    }
)
