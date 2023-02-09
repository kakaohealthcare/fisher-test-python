/*
This version based on the version in http://netlib.org/toms/
643.gz

By sector.rest@kakaohealthcare.com
*/
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))


int f2xact(
    int nrow,
    int ncol,
    const int table[],
    int ldtabl,
    double expect,
    double percnt,
    double emin,
    double* prt,
    double* pre,
    double* fact,
    int* ico,
    int* iro,
    int* kyy,
    int* idif,
    int* irn,
    int* key,
    int ldkey,
    int* ipoin,
    double* stp,
    int ldstp,
    int* ifrq,
    double* dlp,
    double* dsp,
    double* tm,
    int* key2,
    int* iwk,
    double* rwk);
int f3xact(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* dlp,
    int* mm, 
    const double fact[],
    int* ico,
    int* iro,
    int* it,
    int* lb,
    int* nr,
    int* nt,
    int* nu,
    int* itc,
    int* ist,
    double* stv,
    double* alen,
    double tol);
void f4xact(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* dsp,
    const double fact[],
    int* icstk,
    int* ncstk,
    int* lstk,
    int* mstk,
    int* nstk,
    int* nrstk,
    int* irstk,
    double* ystk,
    double tol);
int f5xact(
    double pastp,
    double tol,
    int* kval,
    int* key,
    int ldkey,
    int* ipoin,
    double* stp,
    int ldstp,
    int* ifrq,
    int* npoin,
    int* nr,
    int* nl,
    int ifreq,
    int* itop,
    int ipsh);
void f6xact(
    int nrow,
    int* irow,
    int* iflag,
    int* kyy,
    int* key,
    int ldkey,
    int* last,
    int* ipn);
void f7xact(
    int nrow,
    int* imax,
    int* idif,
    int* k,
    int* ks,
    int* iflag);
void f8xact(const int irow[], int is, int i1, int izero, int* new);
double f9xact(int n, int mm, const int ir[], const double fact[]);
void f10act(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* val,
    int* xmin,
    const double fact[],
    int* nd,
    int* ne,
    int* m);

void f11act(const int irow[], int i1, int i2, int* new);

int prterr(int code, const char* message);
int iwork(int iwkmax, int* iwkpt, int number, int itype);
void isort(int n, int* ix);
double gammds(double y, double p, int* ifault);
double alogam(double x, int* ifault);

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

  int i, i1, i10, i2, i3, i3a, i3b, i3c, i4, i5, i6, i7, \
      i8, i9, i9a, iiwk, ireal, irwk, iwkmax, iwkpt, \
      j, k, kk, ldkey, ldstp, nco, nro, \
      ntot, numb;

  /* Original comment was
    "
      Set IREAL = 4 for DOUBLE PRECISION
      et IREAL = 3 for SINGLE PRECISION
    "

    but, will fix IREAL to 4
  */
  ireal = 4;

  // iwkmax(workspace) must be even
  iwkmax = 2 * (int)(workspace / 2);

  // create iwkmax * 4 = iwkmax / 2 * 8 bytes size workspace
  double* equivalence = (double*)malloc(iwkmax / 2 * sizeof(double));

#define dwrk(i) equivalence + i - 1
#define iwrk(i) ((int*)equivalence) + i - 1
#define rwrk(i) ((float*)equivalence) + i - 1

  iwkpt = 1;

  if (nrow > ldtabl) {
    return prterr(1, "NROW must be less than or equal to LDTABL.");
  }

  ntot = 0;
  for (i = 0; i < nrow; ++i) {
    for (j = 0; j < ncol; ++j) {
      if (table[j + i * ldtabl] < 0) {
        return prterr(2, "All elements of TABLE must be positive.");
      }
      ntot = ntot + table[j + i * ldtabl];
    }
  }
  if (ntot == 0) {
    return prterr(3, "All elements of TABLE are zero.");
  }

  nco = max(nrow, ncol);
  nro = nrow + ncol - nco;
  k = nrow + ncol + 1;
  kk = k * nco;

#define call_iwork(var, iwkmax, iwkpt, number, type) \
  var = iwork(iwkmax, &iwkpt, number, type); \
  if (var < 0) { return prterr(40, "Out of workspace."); }

  call_iwork(i1, iwkmax, iwkpt, ntot + 1, ireal)
  call_iwork(i2, iwkmax, iwkpt, nco, 2)
  call_iwork(i3, iwkmax, iwkpt, nco, 2)
  call_iwork(i3a, iwkmax, iwkpt, nco, 2)
  call_iwork(i3b, iwkmax, iwkpt, nro, 2)
  call_iwork(i3c, iwkmax, iwkpt, nro, 2)
  call_iwork(iiwk, iwkmax, iwkpt, max(5 * k + 2 * kk, 800 + 7 * max(nrow, ncol)), 2)
  call_iwork(irwk, iwkmax, iwkpt, max(400 + max(nrow, ncol) + 1, k), ireal)

  numb = 18 + 10 * mult;
  ldkey = (iwkmax - iwkpt + 1) / numb;

  ldstp = mult * ldkey;

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
  // Do this for index to start from 1 (to match code with Fortran)
  // fact is excluded, fact index starts from 0 (in Fortran too)
  table -= ldtabl + 1;  // go back 1 row
  --ico; --iro; --kyy; --idif; --irn; --key;
  --ipoin; --stp;
  --ifrq; --dlp; --dsp; --tm; --key2; --iwk; --rwk;

  int i, i31, i310, i311, i32, i33, i34, i35, i36, i37, \
      i38, i39, i41, i42, i43, i44, i45, i46, i47, i48, \
      iflag, ifreq, ii, ikkey, ikstp, ikstp2, ipn, ipo, \
      itmp, itop, itp, j, jkey, jstp, jstp2, jstp3, jstp4, \
      k, k1, kb, kd, kmax, ks, kval, last, n, ncell, nco, \
      nrb, nro, nro2, ntot, ifault, imax;
  double dd, ddf, df, drn, dro, dspt, emn, obs, obs2, \
         obs3, pastp, pv, tmp, tol;
  int chisq, ipsh;

    /*
    IMAX is the largest representable
    integer on the machine
  */
  imax = INT_MAX;
  /*
    TOL is chosen as the square root of
    the smallest relative spacing
  */
  tol = 3.45254e-7;
  /*
    EMX is a large positive value used 
    in comparing expected values
  */
  const double emx = 1.0e30;

  // Initialize KEY array
  for (i = 1; i <= 2 * ldkey; ++i) {
    key[i] = -9999;
    key2[i] = -9999;
  }

  // Initialize parameters
  *pre = 0.0;
  itop = 0;
  if (expect > 0.0) {
    emn = emin;
  } else {
    emn = emx;
  }

  // Initialize pointers for workspace
  k = max(nrow, ncol);
  // f3xact
  i31 = 1;
  i32 = i31 + k;
  i33 = i32 + k;
  i34 = i33 + k;
  i35 = i34 + k;
  i36 = i35 + k;
  i37 = i36 + k;
  i38 = i37 + k;
  i39 = i38 + 400;
  i310 = 1;
  i311 = 401;

  // f4axt
  k = nrow + ncol + 1;
  i41 = 1;
  i42 = i41 + k;
  i43 = i42 + k;
  i44 = i43 + k;
  i45 = i44 + k;
  i46 = i45 + k;
  i47 = i46 + k * max(nrow, ncol);
  i48 = 1;

  // Check table dimensions
  if (nrow > ldtabl) {
    return prterr(1, "NROW must be less than or equal to LDTABL.");
  }
  if (ncol <= 1) {
    return prterr(4, "NCOL must be greater than 1.0.");
  }

#define get_table(i, j) table[(j) + (i) * ldtabl]

  // Compute row marginals and total
  ntot = 0;
  for (i = 1; i <= nrow; ++i) {
    iro[i] = 0;
    for (j = 1; j <= ncol; ++j) {
      if (get_table(i, j) < -0.0001) {
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
  isort(nrow, &iro[1]);
  isort(ncol, &ico[1]);

  // Determine row and column marginals
  if (nrow > ncol) {
    nro = ncol;
    nco = nrow;
    // Interchange row and column marginals
    for (i = 1; i <= nrow; ++i) {
      itmp = iro[i];
      if (i <= ncol) iro[i] = ico[i];
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
    kmax = (iro[nro] + 1) * kyy[nro - 1];
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
    fact[i] = fact[i - 1] + log((double)i);
    j = i + 1;
    if (j <= ntot) {
      fact[j] = fact[i] + fact[2] + fact[j / 2] - fact[j / 2 - 1];
    }
  }

  // Compute observed path length: OBS
  obs = tol;
  ntot = 0;
  for (j = 1; j <= nco; ++j) {
    dd = 0.0;
    for (i = 1; i <= nro; ++i) {
      if (nrow <= ncol) {
        dd += fact[get_table(i, j)];
        ntot += get_table(i, j);
      } else {
        dd += fact[get_table(j, i)];
        ntot += get_table(j, i);
      }
    }
    obs = obs + fact[ico[j]] - dd;
  }

  // Denominator of observed table: DRO
  dro = f9xact(nro, ntot, &iro[1], fact);
  *prt = exp(obs - dro);

  // Initialize pointers
  k = nco;
  last = ldkey + 1;
  jkey = ldkey + 1;
  jstp = ldstp + 1;
  jstp2 = 3 * ldstp + 1;
  jstp3 = 4 * ldstp + 1;
  jstp4 = 5 * ldstp + 1;
  ikkey = 0;
  ikstp = 0;
  ikstp2 = 2 * ldstp;
  ipo = 1;
  ipoin[1] = 1;
  stp[1] = 0.0;
  ifrq[1] = 1;
  ifrq[ikstp2 + 1] = -1;

L110:
  kb = nco - k + 1;
  ks = 0;
  n = ico[kb];
  kd = nro + 1;
  kmax = nro;

  // IDIF is the difference in going to th daughter
  for (i = 1; i <= nro; ++i) {
    idif[i] = 0;
  }

L130:
  // Generate the first daughter
  --kd;
  ntot = min(n, iro[kd]);
  idif[kd] = ntot;
  if (idif[kmax] == 0) --kmax;
  n -= ntot;
  if (n > 0 && kd != 1) goto L130;
  if (n != 0) goto L310;

  k1 = k - 1;
  n = ico[kb];
  ntot = 0;
  for (i = kb + 1; i <= nco; ++i) {
    ntot += ico[i];
  }

L150:
  // Arc to daughter length=ICO(KB)
  for (i = 1; i <= nro; ++i) {
    irn[i] = iro[i] - idif[i];
  }

  // Sort irn
  if (k1 > 1) {
    if (nro == 2) {
      if (irn[1] > irn[2]) {
        ii = irn[1];
        irn[1] = irn[2];
        irn[2] = ii;
      }
    } else if (nro == 3) {
      ii = irn[1];
      if (ii > irn[3]) {
        if (ii > irn[2]) {
          if (irn[2] > irn[3]) {
            irn[1] = irn[3];
            irn[3] = ii;
          } else {
            irn[1] = irn[2];
            irn[2] = irn[3];
            irn[3] = ii;
          }
        } else {
          irn[1] = irn[3];
          irn[3] = irn[2];
          irn[2] = ii;
        }
      } else if (ii > irn[2]) {
        irn[1] = irn[2];
        irn[2] = ii;
      } else if (irn[2] > irn[3]) {
        ii = irn[2];
        irn[2] = irn[3];
        irn[3] = ii;
      }
    } else {
      for (j = 2; j <= nro; ++j) {
        i = j - 1;
        ii = irn[j];
L170:
        if (ii < irn[i]) {
          irn[i + 1] = irn[i];
          --i;
          if (i > 0) goto L170;
        }
        irn[i + 1] = ii;
      }
    }
    // Adjust start for zero
    for (i = 1; i <= nro; ++i) {
      if (irn[i] != 0) {
        goto L200;
      }
    }
L200:
    nrb = i;
    nro2 = nro - i + 1;
  } else {
    nrb = 1;
    nro2 = nro;
  }

  // Some table values
  ddf = f9xact(nro, n, &idif[1], fact);
  drn = f9xact(nro2, ntot, &irn[nrb], fact) - dro + ddf;

  // Get hash value
  if (k1 > 1) {
    kval = irn[1] + irn[2] * kyy[2];
    for (i=3; i <= nro; ++i) {
      kval += irn[i] * kyy[i];
    }

    // Get hash table entry
    i = kval % (2 * ldkey) + 1;

    // Search for unused location
    for (itp = i; itp <= 2 * ldkey; ++itp) {
      ii = key2[itp];
      if (ii == kval) goto L240;
      else if (ii < 0) {
        key2[itp] = kval;
        dlp[itp] = 1.0;
        dsp[itp] = 1.0;
        goto L240;
      }
    }

    for (itp = 1; itp <= i - 1; ++itp) {
      ii = key2[itp];
      if (ii == kval) goto L240;
      else if (ii < 0) {
        key2[itp] = kval;
        dlp[itp] = 1.0;
        goto L240;
      }
    }

    return prterr(
        6,
        "LDKEY is too small.  It is not possible to \
        give thevalue of LDKEY required, but you could \
        try doubling LDKEY (and possibly LDSTP)"
    );
  }

L240:
  // boolean value
  ipsh = 1;

  // Recover pastp
  ipn = ipoin[ipo + ikkey];
  pastp = stp[ipn + ikstp];
  ifreq = ifrq[ipn + ikstp];

  // Compute shortest and longest path
  if (k1 > 1) {
    obs2 = obs - fact[ico[kb + 1]] - fact[ico[kb + 2]] - ddf;
    for (i = 3; i <= k1; ++i) {
      obs2 -= fact[ico[kb + i]];
    }

    if (dlp[itp] > 0.0) {
      dspt = obs - obs2 - ddf;
      // Compute longest path
      dlp[itp] = 0.0;

      int f3_ret = f3xact(
          nro2,
          &irn[nrb],
          k1,
          &ico[kb + 1],
          &dlp[itp],
          &ntot,
          fact,
          &iwk[i31],
          &iwk[i32],
          &iwk[i33],
          &iwk[i34],
          &iwk[i35],
          &iwk[i36],
          &iwk[i37],
          &iwk[i38],
          &iwk[i39],
          &rwk[i310],
          &rwk[i311],
          tol
      );
      if (f3_ret != 0) return f3_ret;

      dlp[itp] = min(0.0, dlp[itp]);

      // Compute shortest path
      dsp[itp] = dspt;
      f4xact(
          nro2,
          &irn[nrb],
          k1,
          &ico[kb + 1],
          &dsp[itp],
          fact,
          &iwk[i47],
          &iwk[i41],
          &iwk[i42],
          &iwk[i43],
          &iwk[i44],
          &iwk[i45],
          &iwk[i46],
          &rwk[i48],
          tol
      );
      dsp[itp] = min(0.0, dsp[itp] - dspt);

      // Use chi-squared approximation?
      if ((double)(irn[nrb] * ico[kb + 1]) / (double)ntot > emn) {
        ncell = 0.0;
        for (i = 1; i <= nro2; ++i) {
          for (j = 1; j <= k1; ++j) {
            if (irn[nrb + i - 1] * ico[kb + j] >= ntot * expect) {
              ++ncell;
            }
          }
        }
        if (ncell * 100 >= k1 * nro2 * percnt) {
          tmp = 0.0;
          for (i = 1; i <= nro2; ++i) {
            tmp += fact[irn[nrb + i - 1]] - fact[irn[nrb + i - 1] - 1];
          }
          tmp *= (k1 - 1);
          for (j= 1; j <= k1; ++j) {
            tmp += (nro2 - 1) * (fact[ico[kb + j]] - fact[ico[kb + j] - 1]);
          }

          df = (nro2 - 1) * (k1 - 1);
          tmp += df * 1.83787706640934548356065947281;
          tmp -= (nro2 * k1 - 1) * (fact[ntot] - fact[ntot-1]);
          tm[itp] = -2.0 * (obs - dro) - tmp;
        } else {
          // tm(itp) set to a flag value
          tm[itp] = -9876.0;
        }
      } else {
        tm[itp] = -9876.0;
      }
    }

    obs3 = obs2 - dlp[itp];
    obs2 = obs2 - dsp[itp];
    if (tm[itp] == -9876.0) {
      chisq = 0;
    } else {
      chisq = 1;
      tmp = tm[itp];
    }
  } else {
    obs2 = obs - drn - dro;
    obs3 = obs2;
  }

L300:
  // Process node with new PASTP
  if (pastp <= obs3) {
    // Update pre
    *pre += (double)ifreq * exp(pastp + drn);
  } else if (pastp < obs2) {
    if (chisq) {
      df = (nro2 - 1) * (k1 - 1);
      pv = 1.0 - gammds(
          max(0.0, tmp + 2.0 * (pastp + drn)) / 2.0,
          df / 2.0,
          &ifault);
      *pre += (double)ifreq *exp(pastp + drn) * pv;
    } else {
      // Put daughter on queue
      int f5_ret = f5xact(
          pastp + ddf,
          tol,
          &kval,
          &key[jkey],
          ldkey,
          &ipoin[jkey],
          &stp[jstp],
          ldstp,
          &ifrq[jstp],
          &ifrq[jstp2],
          &ifrq[jstp3],
          &ifrq[jstp4],
          ifreq,
          &itop,
          ipsh
      );
      if (f5_ret != 0) return f5_ret;
      ipsh = 0;
    }
  }

  // Get next PASTP on chain
  ipn = ifrq[ipn + ikstp2];
  if (ipn > 0) {
    pastp = stp[ipn + ikstp];
    ifreq = ifrq[ipn + ikstp];
    goto L300;
  }

  // Generate a new daughter node
  f7xact(
      kmax,
      &iro[1],
      &idif[1],
      &kd,
      &ks,
      &iflag
  );
  if (iflag != 1) goto L150;

L310:
  // Go get a new mother from stage K
  iflag = 1;
  f6xact(
      nro,
      &iro[1],
      &iflag,
      &kyy[1],
      &key[ikkey + 1],
      ldkey,
      &last,
      &ipo
  );
  // Update pointers
  if (iflag == 3) {
    --k;
    itop = 0;
    ikkey = jkey - 1;
    ikstp = jstp - 1;
    ikstp2 = jstp2 - 1;
    jkey = ldkey - jkey + 2;
    jstp = ldstp - jstp + 2;
    jstp2 = 2 * ldstp + jstp;
    for (i = 1; i <= 2 * ldkey; ++i) {
      key2[i] = -9999;
    }
    if (k >= 2) {
      goto L310;
    }
  } else {
    goto L110;
  }

  return 0;
}

int f3xact(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* dlp,
    int* mm, 
    const double fact[],
    int* ico,
    int* iro,
    int* it,
    int* lb,
    int* nr,
    int* nt,
    int* nu,
    int* itc,
    int* ist,
    double* stv,
    double* alen,
    double tol) {
/*
-----------------------------------------------------------------------
  Name:       F3XACT

  Purpose:    Computes the shortest path length for a given table.

  Usage:      CALL F3XACT (NROW, IROW, NCOL, ICOL, DLP, MM, FACT, ICO,
                          IRO, IT, LB, NR, NT, NU, ITC, IST, STV, ALEN,
                          TOL)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     IROW   - Vector of length NROW containing the row sums for the
              table.  (Input)
     NCOL   - The number of columns in the table.  (Input)
     ICOL   - Vector of length K containing the column sums for the
              table.  (Input)
     DLP    - The longest path for the table.  (Output)
     MM     - The total count in the table.  (Output)
     FACT   - Vector containing the logarithms of factorials.  (Input)
     ICO    - Work vector of length MAX(NROW,NCOL).
     IRO    - Work vector of length MAX(NROW,NCOL).
     IT     - Work vector of length MAX(NROW,NCOL).
     LB     - Work vector of length MAX(NROW,NCOL).
     NR     - Work vector of length MAX(NROW,NCOL).
     NT     - Work vector of length MAX(NROW,NCOL).
     NU     - Work vector of length MAX(NROW,NCOL).
     ITC    - Work vector of length 400.
     IST    - Work vector of length 400.
     STV    - Work vector of length 400.
     ALEN   - Work vector of length MAX(NROW,NCOL).
     TOL    - Tolerance.  (Input)
-----------------------------------------------------------------------
*/
  int i, ic1, ic2, ii, ipn, irl, itp, k, key, ks, kyy, lev;
  int n11, n12, nc1, nc1s, nco, nct, nn, nn1, nr1, nro, nrt;
  double v, val, vmn;
  int xmin;

  const int ldst = 200;
  int nst = 0;
  int nitc = 0;

  // Do this for index to start from 1 (to match code with Fortran)
  // fact, alen are excluded, (starts from 0)
  --irow; --icol;
  --ico; --iro; --it; --lb; --nr; --nt; --nu; --itc; --ist;
  --stv;

  for (i = 0; i <= ncol; ++i) {
    alen[i] = 0.0;
  }
  for (i = 1; i <= 400; ++i) {
    ist[i] = -1;
  }

  // nrow is 1
  if (nrow <= 1) {
    if (nrow > 0) {
      *dlp -= fact[icol[1]];
      for (i = 2; i <= ncol; ++i) {
        *dlp -= fact[icol[i]];
      }
    }
    goto L9000;
  }

  // ncol is 1
  if (ncol <= 1) {
    if (ncol > 0) {
      *dlp -= fact[irow[1]] + fact[irow[2]];
      for (i = 3; i <= nrow; ++i) {
        *dlp -= fact[irow[i]];
      }
    }
    goto L9000;
  }

  // 2 by 2 table
  if (nrow * ncol == 4) {
    n11 = (irow[1] + 1) * (icol[1] + 1) /  (*mm + 2);
    n12 = irow[1] - n11;
    *dlp -= (fact[n11] + fact[n12] + fact[icol[1] - n11] + fact[icol[2]-n12]);
    goto L9000;
  }

  // Test for optimal table
  val = 0.0;
  xmin = 0;
  if (irow[nrow] <= irow[1] + ncol) {
    f10act(
        nrow,
        &irow[1],
        ncol,
        &icol[1],
        &val,
        &xmin,
        fact,
        &lb[1],
        &nu[1],
        &nr[1]
    );
  }
  if (xmin == 0) {
    if (icol[ncol] <= icol[1] + nrow) {
      f10act(
          ncol,
          &icol[1],
          nrow,
          &irow[1],
          &val,
          &xmin,
          fact,
          &lb[1],
          &nu[1],
          &nr[1]
      );
    }
  }

  if (xmin) {
    *dlp -= val;
    goto L9000;
  }

  // Setup for dynamic programming
  nn = *mm;

  // Minimize ncol
  if (nrow >= ncol) {
    nro = nrow;
    nco = ncol;
    for (i = 1; i <= nrow; ++i) {
      iro[i] = irow[i];
    }

    ico[1] = icol[1];
    nt[1] = nn - ico[1];
    for (i = 2; i <= ncol; ++i) {
      ico[i] = icol[i];
      nt[i] = nt[i - 1] - ico[i];
    }
  } else {
    nro = ncol;
    nco = nrow;

    ico[1] = irow[1];
    nt[1] = nn - ico[1];
    for (i = 2; i <= nrow; ++i) {
      ico[i] = irow[i];
      nt[i] = nt[i - 1] - ico[i];
    }

    for (i =1; i <= ncol; ++i) {
      iro[i] = icol[i];
    }
  }

  // Initialize pointers
  vmn = 1.0e10;
  nc1s = nco - 1;
  irl = 1;
  ks = 0;
  k = ldst;
  kyy = ico[nco] + 1;
  goto L100;

  // Test for optimality
L90:
  xmin = 0;
  if (iro[nro] <= iro[irl] + nco) {
    f10act(
        nro,
        &iro[irl],
        nco,
        &ico[1],
        &val,
        &xmin,
        fact,
        &lb[1],
        &nu[1],
        &nr[1]
    );
  }
  if (xmin == 0) {
    if (ico[nco] <= ico[1] + nro) {
      f10act(
          nco,
          &ico[1],
          nro,
          &iro[irl],
          &val,
          &xmin,
          fact,
          &lb[1],
          &nu[1],
          &nr[1]
      );
    }
  }

  if (xmin == 1) {
    if (val < vmn) vmn = val;
    goto L200;
  }

L100:
  // Setup to generate new node
  lev = 1;
  nr1 = nro - 1;
  nrt = iro[irl];
  nct = ico[1];
  lb[1] = (int)((double)((nrt + 1) * (nct + 1)) / (double)(nn + nr1 * nc1s + 1) - tol) - 1;
  nu[1] = (int)((double)((nrt + nc1s) * (nct + nr1)) / (double)(nn + nr1 + nc1s)) - lb[1] + 1;
  nr[1] = nrt - lb[1];

L110:
  // Generate a node
  --nu[lev];
  if (nu[lev] == 0) {
    if (lev == 1) goto L200;
    --lev;
    goto L110;
  }
  ++lb[lev];
  --nr[lev];

L120:
  alen[lev] = alen[lev - 1] + fact[lb[lev]];
  if (lev < nc1s) {
    nn1 = nt[lev];
    nrt = nr[lev];
    ++lev;
    nc1 = nco - lev;
    nct = ico[lev];
    lb[lev] = (double)((nrt + 1) * (nct + 1)) / (double)(nn1 + nr1 * nc1 + 1) - tol;
    nu[lev] = (double)((nrt + nc1) * (nct + nr1)) / (double)(nn1 + nr1 + nc1) - lb[lev] + 1;
    nr[lev] = nrt - lb[lev];
    goto L120;
  }
  alen[nco] = alen[lev] + fact[nr[lev]];
  lb[nco] = nr[lev];

  v = val + alen[nco];
  if (nro == 2) {
    // Only 1 row left;
    v += fact[ico[1] - lb[1]] + fact[ico[2] - lb[2]];
    for (i = 3; i <= nco; ++i) {
      v += fact[ico[i] - lb[i]];
    }
    if (v < vmn) vmn = v;
  } else if (nro == 3 && nco == 2) {
    // 3 rows and 2 columns
    nn1 = nn - iro[irl] + 2;
    ic1 = ico[1] - lb[1];
    ic2 = ico[2] - lb[2];
    n11 = (iro[irl + 1] + 1) * (ic1 + 1) / nn1;
    n12 = iro[irl + 1] - n11;
    v += fact[n11] + fact[n12] + fact[ic1 - n11] + fact[ic2 - n12];
    if (v < vmn) vmn = v;
  } else {
    // Column marginals are new node
    for (i = 1; i <= nco; ++i) {
      it[i] = ico[i] - lb[i];
    }

    // Sort column marginals
    if (nco == 2) {
      if (it[1] > it[2]) {
        ii = it[1];
        it[1] = it[2];
        it[2] = ii;
      }
    }
    else if (nco == 3) {
      ii = it[1];
      if (ii > it[3]) {
        if (ii > it[2]) {
          if (it[2] > it[3])  {
            it[1] = it[3];
            it[3] = ii;
          } else {
            it[1] = it[2];
            it[2] = it[3];
            it[3] = ii;
          }
        } else {
          it[1] = it[3];
          it[3] = it[2];
          it[2] = ii;
        }
      } else if (ii > it[2]) {
        it[1] = it[2];
        it[2] = ii;
      } else if (it[2] > it[3]) {
        ii = it[2];
        it[2] = it[3];
        it[3] = ii;
      }
    } else {
      isort(nco, &it[1]);
    }

    // Compute hash value
    key = it[1] * kyy + it[2];
    for (i = 3; i <= nco; ++i) {
      key = it[i] + key * kyy;
    }

    // Table index
    ipn = key % ldst + 1;

    // Find empty position
    ii = ks + ipn;
    for (itp=ipn; itp <= ldst; ++itp) {
      if (ist[ii] < 0) {
        goto L180;
      } else if (ist[ii] == key) {
        goto L190;
      }
      ++ii;
    }

    ii = ks + 1;
    for (itp = 1; itp <= ipn - 1; ++itp) {
      if (ist[ii] < 0) {
        goto L180;
      } else if (ist[ii] == key) {
        goto L190;
      }
      ++ii;
    }

    return prterr(30, "Stack length exceed in f3xact. This problem should not occur.");

L180:
    // Push onto tack
    ist[ii] = key;
    stv[ii] = v;
    ++nst;
    ii = nst + ks;
    itc[ii] = itp;
    goto L110;

L190:
    // Marginals already on stack
    stv[ii] = min(v, stv[ii]);
  }
  goto L110;

L200:
  // Pop item from stack
  if (nitc > 0) {
    // Stack index
    itp = itc[nitc + k] + k;
    --nitc;
    val = stv[itp];
    key = ist[itp];
    ist[itp] = -1;

    // Compute marginals
    for (i = nco; i >= 2; --i) {
      ico[i] = key % kyy;
      key /= kyy;
    }
    ico[1] = key;
    // Set up nt array
    nt[1] = nn - ico[1];
    for (i = 2; i <= nco; ++i) {
      nt[i] = nt[i - 1] - ico[i];
    }
    goto L90;
  } else if (nro > 2 && nst > 0) {
    // Go to next level
    nitc = nst;
    nst = 0;
    k = ks;
    ks = ldst - ks;
    nn -= iro[irl];
    ++irl;
    --nro;
    goto L200;
  }

  *dlp -= vmn;

L9000:

  return 0;
}

void f4xact(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* dsp,
    const double fact[],
    int* icstk,
    int* ncstk,
    int* lstk,
    int* mstk,
    int* nstk,
    int* nrstk,
    int* irstk,
    double* ystk,
    double tol) {
/*
-----------------------------------------------------------------------
  Name:       F4XACT

  Purpose:    Computes the longest path length for a given table.

  Usage:      CALL F4XACT (NROW, IROW, NCOL, ICOL, DSP, FACT, ICSTK,
                          NCSTK, LSTK, MSTK, NSTK, NRSTK, IRSTK, YSTK,
                          TOL)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     IROW   - Vector of length NROW containing the row sums for the
              table.  (Input)
     NCOL   - The number of columns in the table.  (Input)
     ICOL   - Vector of length K containing the column sums for the
              table.  (Input)
     DSP    - The shortest path for the table.  (Output)
     FACT   - Vector containing the logarithms of factorials.  (Input)
     ICSTK  - NCOL by NROW+NCOL+1 work array.
     NCSTK  - Work vector of length NROW+NCOL+1.
     LSTK   - Work vector of length NROW+NCOL+1.
     MSTK   - Work vector of length NROW+NCOL+1.
     NSTK   - Work vector of length NROW+NCOL+1.
     NRSTK  - Work vector of length NROW+NCOL+1.
     IRSTK  - NROW by MAX(NROW,NCOL) work array.
     YSTK   - Work vector of length NROW+NCOL+1.
     TOL    - Tolerance.  (Input)
-----------------------------------------------------------------------
*/
  int i, ic1, ict, ir1, irt, istk, j, k, l, m, mn, n, nco, nro;
  double amx, y;

  // Do this for index to start from 1 (to match code with Fortran)
  --irow; --icol;
  --ncstk; --lstk; --mstk; --nstk; --nrstk;
  --ystk;
  icstk -= ncol + 1;
  irstk -= nrow + 1;

#define get_irstk(i, j) irstk[i + (j) * nrow]
#define get_icstk(i, j) icstk[i + (j) * ncol]

  // Take care of the easy cases firstkt
  if (nrow == 1) {
    for (i = 1; i <= ncol; ++i) {
      *dsp -= fact[icol[i]];
    }
    goto L9000;
  }

  if (ncol == 1) {
    for (i = 1; i <= nrow; ++i) {
      *dsp -= fact[irow[i]];
    }
    goto L9000;
  }

  if (nrow * ncol == 4) {
    if (irow[2] <= icol[2]) {
      *dsp -= fact[irow[2]] + fact[icol[1]] + fact[icol[2] - irow[2]];
    } else {
      *dsp -= fact[icol[2]] + fact[irow[1]] + fact[irow[2] - icol[2]];
    }
    goto L9000;
  }

  // initialization before loop
  for (i = 1; i <= nrow; ++i) {
    get_irstk(i, 1) = irow[nrow - i + 1];
  }
  for (j = 1; j <= ncol; ++j) {
    get_icstk(j, 1) = icol[ncol - j + 1];
  }

  nro = nrow;
  nco = ncol;
  nrstk[1] = nro;
  ncstk[1] = nco;
  ystk[1] = 0.0;
  y = 0.0;
  istk = 1;
  l = 1;
  amx = 0.0;

L50:
  ir1 = get_irstk(1, istk);
  ic1 = get_icstk(1, istk);
  if (ir1 > ic1) {
    if (nro >= nco) {
      m = nco - 1;
      n = 2;
    } else {
      m = nro;
      n = 1;
    }
  } else if (ir1 < ic1) {
    if (nro <= nco) {
      m = nro - 1;
      n = 1;
    } else {
      m = nco;
      n = 2;
    }
  } else {
    if (nro <= nco) {
      m = nro - 1;
      n = 1;
    } else {
      m = nco - 1;
      n = 2;
    }
  }

L60:
  if (n == 1) {
    i = l;
    j = 1;
  } else {
    i = 1;
    j = l;
  }

  irt = get_irstk(i, istk);
  ict = get_icstk(j, istk);
  mn = irt;
  if (mn > ict) mn = ict;
  y += fact[mn];
  if (irt == ict) {
    --nro;
    --nco;
    f11act(&get_irstk(1, istk), i, nro, &get_irstk(1, istk + 1));
    f11act(&get_icstk(1, istk), j, nco, &get_icstk(1, istk + 1));
  } else if (irt > ict) {
    --nco;
    f11act(&get_icstk(1, istk), j, nco, &get_icstk(1, istk + 1));
    f8xact(&get_irstk(1, istk), irt - ict, i, nro, &get_irstk(1, istk + 1));
  } else {
    --nro;
    f11act(&get_irstk(1, istk), i, nro, &get_irstk(1, istk + 1));
    f8xact(&get_icstk(1, istk), ict - irt, j, nco, &get_icstk(1, istk + 1));
  }

  if (nro == 1) {
    for (k = 1; k <= nco; ++k) {
      y += fact[get_icstk(k, istk + 1)];
    }
    goto L90;
  }

  if (nco == 1) {
    for (k = 1; k <= nro; ++k) {
      y += fact[get_irstk(k, istk + 1)];
    }
    goto L90;
  }

  lstk[istk] = l;
  mstk[istk] = m;
  nstk[istk] = n;
  ++istk;
  nrstk[istk] = nro;
  ncstk[istk] = nco;
  ystk[istk] = y;
  l = 1;
  goto L50;

L90:
  if (y > amx) {
    amx = y;
    if (*dsp - amx <= tol) {
      *dsp = 0.0;
      goto L9000;
    }
  }

L100:
  --istk;
  if (istk == 0) {
    *dsp -= amx;
    if (*dsp - amx <= tol) *dsp = 0.0;
    goto L9000;
  }
  l = lstk[istk] + 1;

L110:
  if (l > mstk[istk]) goto L100;
  n = nstk[istk];
  nro = nrstk[istk];
  nco = ncstk[istk];
  y = ystk[istk];
  if (n == 1) {
    if (get_irstk(l, istk) < get_irstk(l - 1, istk)) goto L60;
  } else if (n == 2) {
    if (get_icstk(l, istk) < get_icstk(l - 1, istk)) goto L60;
  }

  ++l;
  goto L110;

L9000:
  return;
}

int f5xact(
    double pastp,
    double tol,
    int* kval,
    int* key,
    int ldkey,
    int* ipoin,
    double* stp,
    int ldstp,
    int* ifrq,
    int* npoin,
    int* nr,
    int* nl,
    int ifreq,
    int* itop,
    int ipsh) {
/*
-----------------------------------------------------------------------
  Name:       F5XACT

  Purpose:    Put node on stack in network algorithm.

  Usage:      CALL F5XACT (PASTP, TOL, KVAL, KEY, LDKEY, IPOIN, STP,
                          LDSTP, IFRQ, NPOIN, NR, NL, IFREQ, ITOP,
                          IPSH)

  Arguments:
     PASTP  - The past path length.  (Input)
     TOL    - Tolerance for equivalence of past path lengths.  (Input)
     KVAL   - Key value.  (Input)
     KEY    - Vector of length LDKEY containing the key values.
              (Input/output)
     LDKEY  - Length of vector KEY.  (Input)
     IPOIN  - Vector of length LDKEY pointing to the linked list
              of past path lengths.  (Input/output)
     STP    - Vector of length LSDTP containing the linked lists
              of past path lengths.  (Input/output)
     LDSTP  - Length of vector STP.  (Input)
     IFRQ   - Vector of length LDSTP containing the past path
              frequencies.  (Input/output)
     NPOIN  - Vector of length LDSTP containing the pointers to
              the next past path length.  (Input/output)
     NR     - Vector of length LDSTP containing the right object
              pointers in the tree of past path lengths.
              (Input/output)
     NL     - Vector of length LDSTP containing the left object
              pointers in the tree of past path lengths.
              (Input/output)
     IFREQ  - Frequency of the current path length.  (Input)
     ITOP   - Pointer to the top of STP.  (Input)
     IPSH   - Option parameter.  (Input)
              If IPSH is true, the past path length is found in the
              table KEY.  Otherwise the location of the past path
              length is assumed known and to have been found in
              a previous call.
-----------------------------------------------------------------------
*/
  --key; --ipoin; --stp; --ifrq; --npoin; --nr; --nl;

  int ipn, ird, itmp;
  double test1, test2;
  int itp;

  if (ipsh) {
    // Convert KVAL to integer in range 1, ..., LDKEY.
    ird = *kval % ldkey + 1;

    // Search for an unused location
    for (itp = ird; itp <= ldkey; ++itp){
      if (key[itp] == *kval) goto L40;
      if (key[itp] < 0) goto L30;
    }

    for (itp = 1; itp <= ird - 1; ++itp) {
      if (key[itp] == *kval) goto L40;
      if (key[itp] < 0) goto L30;
    }

    // Return if KEY array is full
    return prterr(6, "LDKEY is too small for this problem. \
      It is not possible to estimate the value of LDKEY \
      required, but twice the current value may be sufficient.");

L30:
    // Update KEY
    key[itp] = *kval;
    ++(*itop);
    ipoin[itp] = *itop;

    // Return if STP array full
    if (*itop > ldstp) {
      return prterr(7, "LDSTP is too small for this problem. \
        It is not possible to estimate the value of LDSTP \
        required, but twice the current value may be sufficient");
    }

    // Update STP, etc.
    npoin[*itop] = -1;
    nr[*itop] = -1;
    nl[*itop] = -1;
    stp[*itop] = pastp;
    ifrq[*itop] = ifreq;
    goto L9000;
  }
L40:
  // Find location, if any, of pastp
  ipn = ipoin[itp];
  test1 = pastp - tol;
  test2 = pastp + tol;

L50:
  if (stp[ipn] < test1) {
    ipn = nl[ipn];
    if (ipn > 0) goto L50;
  } else if (stp[ipn] > test2) {
    ipn = nr[ipn];
    if (ipn > 0) goto L50;
  } else {
    ifrq[ipn] += ifreq;
    goto L9000;
  }

  // Return if STP array full
  ++(*itop);
  if (*itop > ldstp) {
    return prterr(7, "LDSTP is too small for this problem. \
        It is not possible to estimate the value of LDSTP \
        required, but twice the current value may be sufficient");
  }

  // Find location to add value
  ipn = ipoin[itp];
  itmp = ipn;

L60:
  if (stp[ipn] < test1) {
    itmp = ipn;
    ipn = nl[ipn];
    if (ipn > 0) {
      goto L60;
    } else {
      nl[itmp] = *itop;
    }
  } else if (stp[ipn] > test2) {
    itmp = ipn;
    ipn = nr[ipn];
    if (ipn > 0) {
      goto L60;
    } else {
      nr[itmp] = *itop;
    }
  }

  // Update STP, etc.
  npoin[*itop] = npoin[itmp];
  npoin[itmp] = *itop;
  stp[*itop] = pastp;
  ifrq[*itop] = ifreq;
  nl[*itop] = -1;
  nr[*itop] = -1;

L9000:

  return 0;
}

void f6xact(
    int nrow,
    int* irow,
    int* iflag,
    int* kyy,
    int* key,
    int ldkey,
    int* last,
    int* ipn) {
/*
-----------------------------------------------------------------------
  Name:       F6XACT

  Purpose:    Pop a node off the stack.

  Usage:      CALL F6XACT (NROW, IROW, IFLAG, KYY, KEY, LDKEY, LAST,
                          IPN)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     IROW   - Vector of length nrow containing the row sums on output.
              (Output)
     IFLAG  - Set to 3 if there are no additional nodes to process.
              (Output)
     KYY    - Constant mutlipliers used in forming the hash table key.
              (Input)
     KEY    - Vector of length LDKEY containing the hash table keys.
              (Input/output)
     LDKEY  - Length of vector KEY.  (Input)
     LAST   - Index of the last key popped off the stack.
              (Input/output)
     IPN    - Pointer to the linked list of past path lengths.
              (Output)
-----------------------------------------------------------------------
*/
  --irow; --kyy; --key;
  int j, kval;

L10:
  ++(*last);
  if (*last <= ldkey) {
    if (key[*last] < 0) goto L10;

    // Get KVAL from the stack
    kval = key[*last];
    key[*last] = -9999;
    for (j = nrow; j >= 2; --j) {
      irow[j] = kval / kyy[j];
      kval -= irow[j] * kyy[j];
    }
    irow[1] = kval;
    *ipn = *last;
  } else {
    *last = 0;
    *iflag = 3;
  }
}

void f7xact(
    int nrow,
    int* imax,
    int* idif,
    int* k,
    int* ks,
    int* iflag) {
/*
-----------------------------------------------------------------------
  Name:       F7XACT

  Purpose:    Generate the new nodes for given marinal totals.

  Usage:      CALL F7XACT (NROW, IMAX, IDIF, K, KS, IFLAG)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     IMAX   - The row marginal totals.  (Input)
     IDIF   - The column counts for the new column.  (Input/output)
     K      - Indicator for the row to decrement.  (Input/output)
     KS     - Indicator for the row to increment.  (Input/output)
     IFLAG  - Status indicator.  (Output)
              If IFLAG is zero, a new table was generated.  For
              IFLAG = 1, no additional tables could be generated.
-----------------------------------------------------------------------
*/
  --imax; --idif;

  int i, k1, m, mm;

  *iflag = 0;

  // Find node which can be incremented, ks
  if (*ks == 0) {
L10:
    ++(*ks);
    if (idif[*ks] == imax[*ks]) goto L10;
  }

  // Find node to decrement (>ks)
  if (idif[*k] > 0 && *k > *ks) {
    --idif[*k];
L30:
    --(*k);
    if (imax[*k] == 0) goto L30;
    m = *k;

L40:
    // Find node to increment (>= ks)
    if (idif[m] >= imax[m]) {
      --m;
      goto L40;
    }
    ++idif[m];

    // Chang eks
    if (m == *ks) {
      if (idif[m] == imax[m]) *ks = *k;
    }
  } else {
L50:
    // Check for finish
    for (k1 = *k + 1; k1 <= nrow; ++k1) {
      if (idif[k1] > 0) goto L70;
    }
    *iflag = 1;
    goto L9000;

L70:
    // Reallocate counts
    mm = 1;
    for (i = 1; i <= *k; ++i) {
      mm += idif[i];
      idif[i] = 0;
    }
    *k = k1;
L90:
    --(*k);
    m = min(mm, imax[*k]);
    idif[*k] = m;
    mm -= m;
    if (mm > 0 && *k != 1) goto L90;

    // Check that all counts reallocated
    if (mm > 0) {
      if (k1 != nrow) {
        *k = k1;
        goto L50;
      }
      *iflag = 1;
      goto L9000;
    }    

    // Get ks
    --idif[k1];
    *ks = 0;
L100:
    ++(*ks);
    if (*ks > *k) goto L9000;
    if (idif[*ks] >= imax[*ks]) goto L100;
  }

L9000:
  return;
}

void f8xact(const int irow[], int is, int i1, int izero, int* new) {
/*
-----------------------------------------------------------------------
  Name:       F8XACT

  Purpose:    Routine for reducing a vector when there is a zero
              element.

  Usage:      CALL F8XACT (IROW, IS, I1, IZERO, NEW)

  Arguments:
     IROW   - Vector containing the row counts.  (Input)
     IS     - Indicator.  (Input)
     I1     - Indicator.  (Input)
     IZERO  - Position of the zero.  (Input)
     NEW    - Vector of new row counts.  (Output)
-----------------------------------------------------------------------
*/
  --irow; --new;
  int i;
  for (i = 1; i <= i1 - 1; ++i) {
    new[i] = irow[i];
  }
  for (i = i1; i <= izero - 1; ++i) {
    if (is >= irow[i + 1]) goto L30;
    new[i] = irow[i + 1];
  }

  i = izero;

L30:
  new[i] = is;

L40:
  ++i;
  if (i > izero) return;
  new[i] = irow[i];
  goto L40;
}

double f9xact(int n, int mm, const int ir[], const double fact[]) {
/*
-----------------------------------------------------------------------
  Name:       F9XACT

  Purpose:    Computes the log of a multinomial coefficient.

  Usage:      F9XACT(N, MM, IR, FACT)

  Arguments:
     N      - Length of IR.  (Input)
     MM     - Number for factorial in numerator.  (Input)
     IR     - Vector of length N containing the numebers for the
              denominator of the factorial.  (Input)
     FACT   - Table of log factorials.  (Input)
     F9XACT  - The log of the multinomal coefficient.  (Output)
-----------------------------------------------------------------------
*/
  --ir;
  int k;
  double ret = fact[mm];
  for (k = 1; k <= n; ++k) {
    ret -= fact[ir[k]];
  }
  return ret;
}

void f10act(
    int nrow,
    const int irow[],
    int ncol,
    const int icol[],
    double* val,
    int* xmin,
    const double fact[],
    int* nd,
    int* ne,
    int* m) {
/*
-----------------------------------------------------------------------
  Name:       F10ACT

  Purpose:    Computes the shortest path length for special tables.

  Usage:      CALL F10ACT (NROW, IROW, NCOL, ICOL, VAL, XMIN, FACT, ND,
                          NE, M)

  Arguments:
     NROW   - The number of rows in the table.  (Input)
     IROW   - Vector of length NROW containing the row totals.  (Input)
     NCOL   - The number of columns in the table.  (Input)
     ICO    - Vector of length NCOL containing the column totals.
              (Input)
     VAL    - The shortest path.  (Output)
     XMIN   - Set to true if shortest path obtained.  (Output)
     FACT   - Vector containing the logarithms of factorials.
              (Input)
     ND     - Workspace vector of length NROW.
     NE     - Workspace vector of length NCOL.
     M      - Workspace vector of length NCOL.

  Chapter:    STAT/LIBRARY Categorical and Discrete Data Analysis
-----------------------------------------------------------------------
*/
  int i, is, ix, nrw1;

  // Do this for index to start from 1 (to match code with Fortran)
  --irow; --icol;
  --nd; --ne; --m;

  for (i = 1; i <= nrow - 1; ++i) {
    nd[i] = 0;
  }

  is = icol[1] / nrow;
  ne[1] = is;
  ix = icol[1] - nrow * is;
  m[1] = ix;
  if (ix != 0) ++nd[ix];

  for (i = 2; i <= ncol; ++i) {
    ix = icol[i] / nrow;
    ne[i] = ix;
    is += ix;
    ix = icol[i] - nrow * ix;
    m[i] = ix;
    if (ix != 0) ++nd[ix];
  }

  for (i = nrow - 2; i >= 1; --i) {
    nd[i] += nd[i+1];
  }

  ix = 0;
  nrw1 = nrow + 1;
  for (i = nrow; i >= 2; --i) {
    ix += is + nd[nrw1 - i] - irow[i];
    if (ix < 0) return;
  }

  for (i = 1; i <= ncol; ++i) {
    ix = ne[i];
    is = m[i];
    *val += is * fact[ix + 1] + (nrow - is) * fact[ix];
  }
  *xmin = true;
}

void f11act(const int irow[], int i1, int i2, int* new) {
/*
-----------------------------------------------------------------------
  Name:       F11ACT

  Purpose:    Routine for revising row totals.

  Usage:      CALL F11ACT (IROW, I1, I2, NEW)

  Arguments:
     IROW   - Vector containing the row totals.  (Input)
     I1     - Indicator.  (Input)
     I2     - Indicator.  (Input)
     NEW    - Vector containing the row totals.  (Input)
-----------------------------------------------------------------------
*/
  --irow; --new;
  int i;
  for (i = 1; i <= i1 - 1; ++i) {
    new[i] = irow[i];
  }

  for (i = i1; i <= i2; ++i) {
    new[i] = irow[i + 1];
  }
}

int prterr(int code, const char *message) {
  printf("%d %s\n", code, message);
  return code;
}

int iwork(int iwkmax, int* iwkpt, int number, int itype) {
/*
c-----------------------------------------------------------------------
  Name:       IWORK

  Purpose:    Routine for allocating workspace.

  Usage:      IWORK (IWKMAX, IWKPT, NUMBER, ITYPE)

  Arguments:
     IWKMAX - Maximum length of workspace.  (Input)
     IWKPT  - Amount of workspace currently allocated.  (Input/output)
     NUMBER - Number of elements of workspace desired.  (Input)
     ITYPE  - Worspace type.  (Input)
              ITYPE  TYPE
                2    Integer
                3    Real
                4    Double Precision
     IWORK  - Index in RWRK, DWRK, or IWRK of the beginning of the
              first element in the workspace array.  (Output)
-----------------------------------------------------------------------
*/
  int ret = *iwkpt;
  if (itype == 2 || itype == 3) {
    *iwkpt += number;
  } else {
    if (ret % 2 != 0) ++ret;
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
  if (i >= j) goto L40;
  kl = i;
  ku = j;
  ikey = i;
  ++j;

L20:
  // Find element in first half
  ++i;
  if (i < j) {
    if (ix[ikey] >=ix[i]) goto L20;
  }

L30:
  // Find element in second half
  --j;
  if (ix[j] > ix[ikey]) goto L30;

  // Interchange
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

double gammds(double y, double p, int* ifault) {
/*
-----------------------------------------------------------------------
  Name:       GAMMDS

  Purpose:    Cumulative distribution for the gamma distribution.

  Usage:      PGAMMA (Q, ALPHA,IFAULT)

  Arguments:
     Q      - Value at which the distribution is desired.  (Input)
     ALPHA  - Parameter in the gamma distribution.  (Input)
     IFAULT - Error indicator.  (Output)
               IFAULT  DEFINITION
                 0     No error
                 1     An argument is misspecified.
                 2     A numerical error has occurred.
     PGAMMA - The cdf for the gamma distribution with parameter alpha
              evaluated at Q.  (Output)
-----------------------------------------------------------------------

       Algorithm AS 147 APPL. Statist. (1980) VOL. 29, P. 113

       Computes the incomplete gamma integral for positive
       parameters Y, P using and infinite series.
*/
  int ifail;
  double a, c, f;
  double ret;
  
  // Checks for the admissibility of arguments and value of F
  *ifault = 1;
  ret = 0.0;
  if (y <= 0.0 || p <= 0.0) return ret;

  *ifault = 2;
  
  /*
    ALOGAM is natural log of gamma function
    no need to test ifail as an error is impossible
  */
  f = exp(p * log(y) - alogam(p + 1.0, &ifail) - y);
  if (f == 0.0) return ret;
  *ifault = 0;

  // Series begins
  c = 1.0;
  ret = 1.0;
  a = p;

L10:
  a += 1.0;
  c = c * y / a;
  ret += c;
  if (c / ret > 1.0e-6) goto L10;
  ret *= f;
  return ret;
}

double alogam(double x, int* ifault) {
/*
-----------------------------------------------------------------------
  Name:       ALOGAM

  Purpose:    Value of the log-gamma function.

  Usage:      ALOGAM (X, IFAULT)

  Arguments:
     X      - Value at which the log-gamma function is to be evaluated.
              (Input)
     IFAULT  - Error indicator.  (Output)
               IFAULT  DEFINITION
                 0     No error
                 1     X .LT. 0
     ALGAMA - The value of the log-gamma function at XX.  (Output)
-----------------------------------------------------------------------

        Algorithm ACM 291, Comm. ACM. (1966) Vol. 9, P. 684

        Evaluates natural logarithm of gamma(x)
        for X greater than zero.

*/
  double a1 = 0.918938533204673;
  double a2 = 0.000595238095238;
  double a3 = 0.000793650793651;
  double a4 = 0.002777777777778;
  double a5 = 0.083333333333333;

  double f, y, z;
  double ret = 0.0;

  *ifault = 1;
  if (x < 0.0) return ret;
  *ifault = 0;

  y = x;
  f = 0.0;
  if (y >= 7.0) goto L30;
  f = y;
L10:
  y += 1.0;
  if (y >= 7.0) goto L20;
  f = f * y;
  goto L10;

L20:
  f = -log(f);
L30:
  z = 1.0 / (y * y);
  ret = f + (y - 0.5) * log(y) - y + a1 + (((-a2 * z + a3) * z - a4) * z + a5) / y;
  return ret;
}
