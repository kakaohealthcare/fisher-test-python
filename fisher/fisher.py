"""
execution code for
quick test
"""

import numpy as np
import ctypes
import os

from ctypes import c_int32, c_long, c_double, byref

lib_path = os.path.join(os.getcwd(), "fisher", "src", "fisher.so")

fisher = ctypes.cdll.LoadLibrary(lib_path)

ary = np.array([
    [1, 3, 1, 4],
    [3, 5, 1, 3]
])

nrow, ncol = ary.shape

ary = np.array(ary.flat).tolist()
print(ary)
ary_c = (c_int32 * len(ary))(*ary)

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
    ncol,
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

# PRT =0.053346 PRE = 0.793284