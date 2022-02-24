/*
This version based on the version in http://netlib.org/toms/

By jungwoo@linewalks.com
*/
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define max(x, y) x > y ? x : y


int f2xact(
    int nrow,
    int ncol,
    const int table[],
    int ldtabl,
    double expect,
    double percnt,
    double emin,
    double *prt,
    double *pre,
    double *fact,
    int *ico,
    int *iro,
    int *kyy,
    int *idif,
    int *irn,
    int *key,
    int ldkey,
    int *ipoin,
    double *stp,
    int ldstp,
    int *ifrq,
    double *dlp,
    double *dsp,
    double *tm,
    int *key2,
    int *iwk,
    double *rwk);

double f9xact(int n, int ntot, const int ir[], const double fact[]);

int prterr(int code, const char *message);
int iwork(int iwkmax, int* iwkpt, int number, int itype);
void isort(int n, int *ix);

int fexact(
    int nrow,
    int ncol,
    const int table[],
    int ldtabl,
    double expect,
    double percnt,
    double emin,
    double *prt,
    double *pre,
    int workspace,
    int mult) {
/*
-----------------------------------------------------------------------
  Name:       FEXACT

  Purpose:    Computes Fisher's exact test probabilities and a hybrid
              approximation to Fisher exact test probabilities for a
              contingency table using the network algorithm.

  Usage:      CALL FEXACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
                          EMIN, PRT, PRE)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     NCOL   - The number of columns in the table.  (Input)
     TABLE  - NROW by NCOL matrix containing the contingency table.
              (Input)
     LDTABL - Leading dimension of TABLE exactly as specified in the
              dimension statement in the calling program.  (Input)
     EXPECT - Expected value used in the hybrid algorithm for
              deciding when to use asymptotic theory probabilities.
              (Input)
              If EXPECT .LE. 0.0 then asymptotic theory probabilities
              are not used and Fisher exact test probabilities are
              computed.  Otherwise, if PERCNT or more of the cells in
              the remaining table have estimated expected values of
              EXPECT or more, with no remaining cell having expected
              value less than EMIN, then asymptotic chi-squared
              probabilities are used.  See the algorithm section of the
              manual document for details.  Use EXPECT = 5.0 to obtain
              the 'Cochran' condition.
     PERCNT - Percentage of remaining cells that must have estimated
              expected  values greater than EXPECT before asymptotic
              probabilities can be used.  (Input)
              See argument EXPECT for details.  Use PERCNT = 80.0 to
              obtain the 'Cochran' condition.
     EMIN   - Minimum cell estimated expected value allowed for
              asymptotic chi-squared probabilities to be used.  (Input)
              See argument EXPECT for details.  Use EMIN = 1.0 to
              obtain the 'Cochran' condition.
     PRT    - Probability of the observed table for fixed marginal
              totals.  (Output)
     PRE    - Table p-value.  (Output)
              PRE is the probability of a more extreme table, where
              'extreme' is in a probabilistic sense.
              If EXPECT .LT. 0 then the Fisher exact probability
              is returned.  Otherwise, an approximation to the
              Fisher exact probability is computed based upon
              asymptotic chi-squared probabilities for ``large''
              table expected values.  The user defines ``large''
              through the arguments EXPECT, PERCNT, and EMIN.

  Remarks:
  1. For many problems one megabyte or more of workspace can be
     required.  If the environment supports it, the user should begin
     by increasing the workspace used to 200,000 units.

  2. In FEXACT, LDSTP = 30*LDKEY.  The proportion of table space used
     by STP may be changed by changing the line MULT = 30 below to
     another value.

  3. FEXACT may be converted to single precision by setting IREAL = 3,
     and converting all DOUBLE PRECISION specifications (except the
     specifications for RWRK, IWRK, and DWRK) to REAL.  This will
     require changing the names and specifications of the intrinsic
     functions ALOG, AMAX1, AMIN1, EXP, and REAL.  In addition, the
     machine specific constants will need to be changed, and the name
     DWRK will need to be changed to RWRK in the call to F2XACT.

  4. Machine specific constants are specified and documented in F2XACT.
     A missing value code is specified in both FEXACT and F2XACT.

  5. Although not a restriction, is is not generally practical to call
     this routine with large tables which are not sparse and in
     which the 'hybrid' algorithm has little effect.  For example,
     although it is feasible to compute exact probabilities for the
     table
            1 8 5 4 4 2 2
            5 3 3 4 3 1 0
           10 1 4 0 0 0 0,
     computing exact probabilities for a similar table which has been
     enlarged by the addition of an extra row (or column) may not be
     feasible.
-----------------------------------------------------------------------
*/

  /*
    AMISS is a missing value indicator
    which is returned when the
    probability is not defined.
  */
  const double amiss = -12345.;

  /* Original comment was
    "
      Set IREAL = 4 for DOUBLE PRECISION
      et IREAL = 3 for SINGLE PRECISION
    "

    but, will fix IREAL to 4
  */
  const int ireal = 4;

  // iwkmax(workspace) must be even
  int iwkmax = 2 * (int)(workspace / 2);

  // create iwkmax * 4 = iwkmax / 2 * 8 bytes size workspace
  double* equivalence = (double*)malloc(iwkmax / 2 * sizeof(double));
  printf("equivalence %p\n", equivalence);

#define dwrk(i) equivalence + i - 1
#define iwrk(i) ((int*)equivalence) + i - 1
#define rwrk(i) ((float*)equivalence) + i - 1

  int iwkpt = 1;

  if (nrow > ldtabl) {
    return prterr(1, "NROW must be less than or equal to LDTABL.");
  }

  int i, j;
  int ntot = 0;

  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      if (table[i + j * ldtabl] < 0) {
        return prterr(2, "All elements of TABLE must be positive.");
      }
      ntot = ntot + table[i + j * ldtabl];
    }
  }
  if (ntot == 0) {
    return prterr(3, "All elements of TABLE are zero.");
  }

  int nco = max(nrow, ncol);
  int nro = nrow + ncol - nco;
  int k = nrow + ncol + 1;
  int kk = k * nco;

#define call_iwork(var, iwkmax, iwkpt, number, type) \
  var = iwork(iwkmax, &iwkpt, number, type); \
  if (var < 0) { return prterr(40, "Out of workspace."); }

  int i1, i2, i3, i3a, i3b, i3c, iiwk, irwk;

  call_iwork(i1, iwkmax, iwkpt, ntot + 1, ireal)
  printf("%d iwkpt %d\n", i1, iwkpt);
  call_iwork(i2, iwkmax, iwkpt, nco, 2)
  printf("%d iwkpt %d\n", i2, iwkpt);
  call_iwork(i3, iwkmax, iwkpt, nco, 2)
  printf("%d iwkpt %d\n", i3, iwkpt);
  call_iwork(i3a, iwkmax, iwkpt, nco, 2)
  call_iwork(i3b, iwkmax, iwkpt, nro, 2)
  call_iwork(i3c, iwkmax, iwkpt, nro, 2)
  call_iwork(iiwk, iwkmax, iwkpt, max(5 * k + 2 * kk, 800 + 7 * max(nrow, ncol)), 2)
  call_iwork(irwk, iwkmax, iwkpt, max(400 + max(nrow, ncol) + 1, k), ireal)

  int numb = 18 + 10 * mult;
  int ldkey = (iwkmax - iwkpt + 1) / numb;

  int ldstp = mult * ldkey;

  int i4, i5, i6, i7, i8, i9, i9a, i10;

  call_iwork(i4, iwkmax, iwkpt, 2 * ldkey, 2)
  call_iwork(i5, iwkmax, iwkpt, 2 * ldkey, 2)
  call_iwork(i6, iwkmax, iwkpt, 2 * ldstp, ireal)
  call_iwork(i7, iwkmax, iwkpt, 6 * ldstp, 2)
  call_iwork(i8, iwkmax, iwkpt, 2 * ldkey, ireal)
  call_iwork(i9, iwkmax, iwkpt, 2 * ldkey, ireal)
  call_iwork(i9a, iwkmax, iwkpt, 2 * ldkey, ireal)
  call_iwork(i10, iwkmax, iwkpt, 2 * ldkey, 2)

  int ret = f2xact(
    nrow,
    ncol,
    table,
    ldtabl,
    expect,
    percnt,
    emin,
    prt,
    pre,
    dwrk(i1),
    iwrk(i2),
    iwrk(i3),
    iwrk(i3a),
    iwrk(i3b),
    iwrk(i3c),
    iwrk(i4),
    ldkey,
    iwrk(i5),
    dwrk(i6),
    ldstp,
    iwrk(i7),
    dwrk(i8),
    dwrk(i9),
    dwrk(i9a),
    iwrk(i10),
    iwrk(iiwk),
    dwrk(irwk)
  );
  printf("equivalence %p\n", equivalence);
  free(equivalence);
  printf("free success %d\n", ret);
  return ret;
}


int f2xact(
    int nrow,
    int ncol,
    const int table[],
    int ldtabl,
    double expect,
    double percnt,
    double emin,
    double *prt,
    double *pre,
    double *fact,
    int *ico,
    int *iro,
    int *kyy,
    int *idif,
    int *irn,
    int *key,
    int ldkey,
    int *ipoin,
    double *stp,
    int ldstp,
    int *ifrq,
    double *dlp,
    double *dsp,
    double *tm,
    int *key2,
    int *iwk,
    double *rwk) {
/*
-----------------------------------------------------------------------
  Name:       F2XACT

  Purpose:    Computes Fisher's exact test for a contingency table,
              routine with workspace variables specified.

  Usage:      CALL F2XACT (NROW, NCOL, TABLE, LDTABL, EXPECT, PERCNT,
                          EMIN, PRT, PRE, FACT, ICO, IRO, KYY, IDIF,
                          IRN, KEY, LDKEY, IPOIN, STP, LDSTP, IFRQ,
                          DLP, DSP, TM, KEY2, IWK, RWK)
-----------------------------------------------------------------------
*/
  /*
    IMAX is the largest representable
    integer on the machine
  */
  const int imax = INT_MAX;
  /*
    TOL is chosen as the square root of
    the smallest relative spacing
  */
  const double tol = 3.45254e-7;
  /*
    EMX is a large positive value used 
    in comparing expected values
  */
  const double emx = 1.0e30;

  int i, j;
  double emn;

  // Do this for index to start from 1 (to match code with Fortran)
  // fact is excluded, fact index starts from 0 (in Fortran too)
  table -= ldtabl + 1;  // go back 1 row
  --ico; --iro; --kyy; --idif; --irn; --key;
  --ipoin; --stp;
  --ifrq; --dlp; --dsp; --tm; --key2; --iwk; --rwk;

  // Initialize KEY array
  for (i = 1; i <= ldkey * 1; ++i) {
    key[i] = -9999;
    key2[i] = -9999;
  }

  // Initialize parameters
  *pre = 0.0;
  int itop = 0;
  if (expect > 0.0) {
    emn = emin;
  } else {
    emn = emx;
  }

  // Initialize pointers for workspace
  int k;

  // f3xact
  k = max(nrow, ncol);
  int i31 = 1;
  int i32 = i31 + k;
  int i33 = i32 + k;
  int i34 = i33 + k;
  int i35 = i34 + k;
  int i36 = i35 + k;
  int i37 = i36 + k;
  int i38 = i37 + k;
  int i39 = i38 + 400;
  int i310 = 1;
  int i311 = 401;

  // f4axt
  k = nrow + ncol + 1;
  int i41 = 1;
  int i42 = i41 + k;
  int i43 = i42 + k;
  int i44 = i43 + k;
  int i45 = i44 + k;
  int i46 = i45 + k;
  int i47 = i46 + k * max(nrow, ncol);
  int i48 = 1;

  // Check table dimensions
  if (nrow > ldtabl) {
    return prterr(1, "NROW must be less than or equal to LDTABL.");
  }
  if (ncol <= 1) {
    return prterr(4, "NCOL must be greater than 1.0.");
  }

#define get_table(i, j) table[i + j * ldtabl]

  // Compute row marginals and total
  int ntot = 0;
  for (i = 1; i <= nrow; ++i) {
    iro[i] = 0;
    for (j = 1; j <= ncol; ++j) {
      if (get_table(i, j) < 0.) {
        return prterr(2, "All elements of TABLE must be positive.");
      }
      iro[i] += get_table(i, j);
      ntot += get_table(i, j);
    }
  }

  if (ntot == 0) {
    return prterr(3, "All elements of TABLE are zero.");
  }

  // Column marginals
  for (i = 1; i <= ncol; ++i) {
    ico[i] = 0;
    for (j = 1; j <= nrow; ++j) {
      ico[i] += get_table(j, i);
    }
  }

  // sort
  isort(nrow, iro);
  isort(ncol, ico);

  // Determine row and column marginals
  int nro, nco, itmp;
  if (nrow > ncol) {
    nro = ncol;
    nco = nrow;
    // Interchange row and column marginals
    for (i = 1; i <= nrow; ++i) {
      itmp = iro[i];
      if (i <= ncol) {
        iro[i] = ico[i];
      }
      ico[i] = itmp;
    }
  } else {
    nro = nrow;
    nco = ncol;
  }

  // Get multiplers for stack
  kyy[1] = 1;
  for (i = 2; i <= nro; ++i) {
    // Hash table multipliers
    if (iro[i - 1] + 1 <= imax / kyy[i - 1]) {
      kyy[i] = kyy[i - 1] * (iro[i - 1] + 1);
      j = j / kyy[i - 1];
    } else {
      return prterr(5, "The hash table key cannot be computed \
        because the largest key is larger than the \
        largest representable integer.  The \
        algorithm cannot proceed.");
    }
  }

  // Maximum product
  if (iro[nro - 1] + 1 < imax / kyy[nro - 1]) {
    // TODO not using kmax maybe?
    int kmax = (iro[nro] + 1) * kyy[nro - 1];
  } else {
    return prterr(5, "The hash table key cannot be computed \
        because the largest key is larger than the \
        largest representable integer.  The \
        algorithm cannot proceed.");
  }

  // Compute log factorials
  fact[0] = 0.0;
  fact[1] = 0.0;
  fact[2] = log(2.0);
  for (i = 3; i <= ntot; i += 2) {
    fact[i] = fact[i - 1] * log((double)i);
    j = i + 1;
    if (j <= ntot) {
      fact[j] = fact[i] + fact[2] + fact[j / 2] - fact[j / 2 - 1];
    }
  }

  // Compute observed path length: OBS
  double obs = tol;
  ntot = 0;
  for (j = 1; j <= nco; ++j) {
    double dd = 0.0;
    for (i = 1; i <= nro; ++i) {
      if (nrow <= ncol) {
        dd = dd + fact[get_table(i, j)];
        ntot = ntot + get_table(i, j);
      } else {
        dd = dd + fact[get_table(j, i)];
        ntot = ntot + get_table(j, i);
      }
    }
    obs = obs + fact[ico[j]] - dd;
  }

  // Denominator of observed table: DRO
  // double dro = f9xact(nro, ntot, iro, fact);
  // *prt = exp(obs - dro);

  // Initialize pointers
  // k = nco;
  // int last = ldkey;
  // int jkey = ldkey;
  // int jstp = ldstp;
  // int jstp2 = 3 * ldstp;
  // int jstp3 = 4 * ldstp;
  // int jstp4 = 5 * ldstp;
  // int ikkey = 0;
  // int ikstp = 0;
  // int ikstp2 = 2 * ldstp;
  // int ipo = 1;
  // ipoin[0] = 1;
  // stp[0] = 0.0;
  // ifrq[0] = 1;
  // ifrq[ikstp2] = -1;

  // int kb = nco - k;
  // int ks = 0;
  // int n = ico[kb];
  // / 


  *prt = 1234.5678;
  *pre = 8765.4321;
  return 0;
}

double f9xact(int n, int ntot, const int ir[], const double fact[]) {
  return 123.456;
}

int prterr(int code, const char *message) {
  printf("%d %s\n", code, message);
  return code;
}

int iwork(int iwkmax, int* iwkpt, int number, int itype) {
  int ret = *iwkpt;
  if (itype == 2 || itype == 3) {
    *iwkpt += number;
  } else {
    if (ret % 2 != 0) {
      ret += 1;
    }
    *iwkpt += 2 * number;
    ret /= 2;
  }
  if (*iwkpt > iwkmax + 1) {
    // Out of workspace
    return -1;
  }
  return ret;
}

void isort(int n, int *ix) {
/*
-----------------------------------------------------------------------
  Name:       ISORT

  Purpose:    Shell sort for an integer vector.

  Usage:      CALL ISORT (N, IX)

  Arguments:
     N      - Lenth of vector IX.  (Input)
     IX     - Vector to be sorted.  (Input/output)
-----------------------------------------------------------------------
*/
  int i, ikey, il[10], it, iu[10], j, kl, ku, m;

  // Do this for index to start from 1 (to match code with Fortran)
  --ix;

  m = 1;
  i = 1;
  j = n;

L10:
  if (i >= j) {
    goto L40;
  }
  kl = i;
  ku = j;
  ikey = i;
  ++j;

// Find element in first half
L20:
  ++i;
  if (i < j) {
    if (ix[ikey] >=ix[i]) {
      goto L20;
    }
  }

// Find element in second half
L30:
  --j;
  if (ix[j] > ix[ikey]) {
    goto L30;
  }
//Interchange
  if (i < j) {
    it = ix[i];
    ix[i] = ix[j];
    ix[j] = it;
    goto L20;
  }
  it = ix[ikey];
  ix[ikey] = ix[j];
  ix[j] = it;

// Save upper and lower subscripts of the array yet to be sorted
  if (m < 11) {
    if (j - kl < ku - j) {
      il[m] = j + 1;
      iu[m] = ku;
      i = kl;
      --j;
    } else {
      il[m] = kl;
      iu[m] = j - 1;
      i = j + 1;
      j = ku;
    }
    ++m;
    goto L10;
  } else {
  }

// Use another segment
L40:
  --m;
  if (m == 0) {
    return;
  }

  i = il[m];
  j = iu[m];
  goto L10;
}

void test_func2(const int table[], int ldtabl) {
  printf("%p\n", table);
  table -= ldtabl;
  printf("%p\n", table);
  printf("%d\n", table[0]);
  printf("%d\n", table[1]);
}

void test_func() {
  // int* equiv = (int*)malloc(5 * sizeof(int));
  // int i;
  // for (i = 0; i < 5; ++i) {
  //   free(&equiv[i]);
  // }

  // const int table[] = {101, 202, 303, 404, 505};
  // test_func2(table, 5);

  // int* a = equiv;
  // printf("%p\n", a);
  // printf("%d\n", *a);

  // --a;
  // printf("%p\n", a);
  // printf("%d\n", *a);
  // free(equiv);
}
