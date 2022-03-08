import os

from setuptools import setup
from distutils.command.build_py import build_py


PACKAGE = "fisher"


class AfterBuildPy(build_py):
  def build_c(self, path, filename):
    obj_path = os.path.join(path, PACKAGE, "src", f"{filename}.so")
    c_path = os.path.join(path, PACKAGE, "src", f"{filename}.c")
    ret = os.system(
        f"gcc -o {obj_path} -shared -fPIC {c_path}"
    )
    if ret != 0:
      raise Exception("gcc build failed.")

  def run(self):
    build_py.run(self)

    installed_path = self.build_lib
    self.build_c(installed_path, "fisher")
    self.build_c(installed_path, "asa159")


setup(
    install_requires=[
        "numpy",
        "scipy"
    ],
    test_suite="tests",
    tests_require=["pytest"],
    packages=["fisher"],
    package_data={
      PACKAGE: ["src/*.c", "src/asa159.so", "src/fisher.so"]
    },
    cmdclass={"build_py": AfterBuildPy},
)
