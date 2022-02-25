"""
execution code for
quick test
"""

import numpy as np
import ctypes
import os

from ctypes import c_long, c_double, byref

lib_path = os.path.join(os.getcwd(), "fisher", "src", "fisher.so")

fisher = ctypes.cdll.LoadLibrary(lib_path)

ary = np.array([
    [100, 3, 100, 204, 404],
    [3, 5, 1, 3, 5]
])

nrow, ncol = ary.shape

ary = np.array(ary.flat).tolist()
ary_c = (c_long * len(ary))(*ary)

expect = c_double(-1.)
percent = c_double(100.)
emin = c_double(0.)

prt = c_double()
pre = c_double()

fisher.test_func()


ret = fisher.fexact(
    nrow,
    ncol,
    ary_c,
    nrow,
    expect,
    percent,
    emin,
    byref(prt),
    byref(pre),
    200000,
    30
)

print(ret)
print(prt.value)
print(pre.value)
