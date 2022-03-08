import os
import sys
import site

from setuptools import (
    setup,
    find_packages,
    Extension
)
from setuptools.command.install import install


PACKAGE = "fisher"


class AfterInstall(install):

  def build_c(self, path, filename):
    obj_path = os.path.join(path, PACKAGE, "src", f"{filename}.so")
    c_path = os.path.join(path, PACKAGE, "src", f"{filename}.c")
    ret = os.system(
        f"gcc -o {obj_path} -shared -fPIC {c_path}"
    )
    if ret != 0:
      raise Exception("gcc build failed.")

  def run(self):
    install.run(self)

    if "--user" in sys.argv:
      installed_path = site.getusersitepackages()
    else:
      site_path_list = site.getsitepackages()
      for site_path in site_path_list:
        if (os.path.exists(os.path.join(site_path, PACKAGE))):
          installed_path = site_path

    self.build_c(installed_path, "fisher")
    self.build_c(installed_path, "asa159")


setup(
    install_requires=[
        "numpy",
        "scipy"
    ],
    test_suite="tests",
    tests_require=["pytest"],
    packages=find_packages(exclude="tests"),
    package_data={
      PACKAGE: ["src/*.c"]
    },
    cmdclass={"install": AfterInstall},
)
