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


def fisher_test(
    ary,
    expect=-1,
    percent=100.,
    emin=0.,
    workspace=200000,
    mult=30
):
  nrow, ncol = ary.shape
  ary = np.array(ary.flat).tolist()
  ary_c = (c_int32 * len(ary))(*ary)

  expect = c_double(-1.)
  percent = c_double(100.)
  emin = c_double(0.)

  prt = c_double()
  pre = c_double()
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
      workspace,
      mult
  )
  if ret != 0:
    raise Exception("fisher_test failed")
  return (prt.value, pre.value)
