import numpy as np
import ctypes
import os

from scipy.special import gammaln

from ctypes import c_int, c_int32, c_double, byref


def load_lib(lib_name):
  
  lib_path = os.path.join(
      os.path.split(__file__)[0],
      "src",
      lib_name
  )
  return ctypes.cdll.LoadLibrary(lib_path)


fisher = load_lib("fisher.so")
asa159 = load_lib("asa159.so")


def nparray_to_carray(ary, type=c_int32):
  return (type * len(ary))(*ary)


def fisher_test(
    ary,
    expect=-1,
    percent=100.,
    emin=0.,
    workspace=200000,
    mult=30,
    simulate_pvalue=False,
    replicate=2000,
    seed=1
):
  if seed == 0:
    raise ValueError("seed must be non-zero")

  if simulate_pvalue:
    sum_row = ary.sum(axis=1)
    sum_col = ary.sum(axis=0)
    ary = ary[sum_row > 0, :][:, sum_col > 0]

    nrow, ncol = ary.shape
    sum_row = ary.sum(axis=1)
    sum_col = ary.sum(axis=0)

    statistics = -np.sum(gammaln(ary + 1))
    results = _simulate_pvalue(nrow, ncol, sum_row, sum_col, replicate, seed)
    almost = 1. + 64. * np.finfo(np.float64).eps

    return (1 + np.sum(results <= (statistics / almost))) / (replicate + 1)
  else:
    nrow, ncol = ary.shape
    ary = ary.ravel().tolist()
    ary = nparray_to_carray(ary)

    expect = c_double(expect)
    percent = c_double(percent)
    emin = c_double(emin)

    prt = c_double()
    pre = c_double()
    ret = fisher.fexact(
        nrow,
        ncol,
        ary,
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
    return pre.value


def _simulate_pvalue(nrow, ncol, sum_row, sum_col, replicate, seed=0):
  sum_row = nparray_to_carray(sum_row)
  sum_col = nparray_to_carray(sum_col)

  results = np.zeros(replicate)

  key = c_int(0)
  seed = c_int(seed)
  matrix = nparray_to_carray(np.zeros(nrow * ncol, dtype=np.int32))

  ntotal = asa159.i4vec_sum(ncol, sum_col)
  fact = nparray_to_carray(np.zeros(ntotal + 1, dtype=np.float64), c_double)
  ierror = c_int(0)
  for idx in range(replicate):
    asa159.rcont2(
        nrow,
        ncol,
        sum_row,
        sum_col,
        byref(key),
        byref(seed),
        fact,
        matrix,
        byref(ierror)
    )
    answer = 0.
    for val in matrix:
      answer -= fact[val]

    results[idx] = answer
  return results
