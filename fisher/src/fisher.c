/*
This version based on the version in http://netlib.org/toms/

By jungwoo@linewalks.com
*/
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

int ptrerr(int code, const char *message);
int iwork(int iwkmax, int* iwkpt, int number, int itype);

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

#define dwrk(i) equivalence + i
#define iwrk(i) ((int*)equivalence) + i
#define rwrk(i) ((float*)equivalence) + i

  // index starts from 0
  int iwkpt = 0;

  if (nrow > ldtabl) {
    return ptrerr(1, "NROW must be less than or equal to LDTABL.");
  }

  int i, j;
  int ntot = 0;

  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      if (table[i + j * ldtabl] < 0) {
        return ptrerr(2, "All elements of TABLE must be positive.");
      }
      ntot = ntot + table[i, j * ldtabl];
    }
  }
  if (ntot == 0) {
    return ptrerr(3, "All elements of TABLE are zero.");
  }

  int nco = max(nrow, ncol);
  int nro = nrow + ncol - nco;
  int k = nrow + ncol + 1;
  int kk = k * nco;

#define call_iwork(var, iwkmax, iwkpt, number, type) \
  var = iwork(iwkmax, &iwkpt, number, type); \
  if (var < 0) { return ptrerr(40, "Out of workspace."); }

  int i1, i2, i3, i3a, i3b, i3c, iiwk, irwk;

  call_iwork(i1, iwkmax, iwkpt, ntot + 1, ireal)
  call_iwork(i2, iwkmax, iwkpt, nco, 2)
  call_iwork(i3, iwkmax, iwkpt, nco, 2)
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
  free(equivalence);
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
  *prt = 1234.5678;
  *pre = 8765.4321;
  return 0;
}


int ptrerr(int code, const char *message) {
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
