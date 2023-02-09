# fisher test python

Python version of Fisher's exact test.

Based on old FORTRAN codes from

- http://netlib.org/toms/
- https://people.sc.fsu.edu/~jburkardt/c_src/asa159/asa159.html

#### Installation

```
pip install fisher-test-python
```

#### Usage

```python
import nump as np
from fisher import fisher_test

ary = np.array([[1,2,3], [2,3,4]])
p_value = fisher_test(ary)

```

## Contact

- Chae, Jungwoo - sector.rest@kakaohealthcare.com
