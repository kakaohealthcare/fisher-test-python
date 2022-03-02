import numpy as np
import pytest
import random

from fisher import fisher_test

from rpy2 import robjects
from rpy2.robjects.numpy2ri import numpy2rpy
from rpy2.robjects.packages import importr


class TestRpy2:
  @pytest.fixture
  def r_func(self):
    stats = importr("stats")
    return stats.fisher_test

  def _run_r_func(self, r_func, ary):
    ary = numpy2rpy(ary)
    return r_func(ary)

  @pytest.mark.parametrize(
      "ary",
      [
          np.random.randint(0, 100, size=(2,2))
          for _ in range(10)
      ]
  )
  def test_fisher_2x2(self, ary, r_func):
    _, pre = fisher_test(ary)
    
    res = self._run_r_func(r_func, ary)
    assert np.isclose(pre, res[0][0])

  @pytest.mark.parametrize(
      "ary",
      [
          np.random.randint(80, 100, size=(2, random.randint(3, 5)))
          for _ in range(10)
      ]
  )
  def test_fisher_2xN(self, ary, r_func):
    _, pre = fisher_test(ary)

    res = self._run_r_func(r_func, ary)
    assert np.isclose(pre, res[0][0])