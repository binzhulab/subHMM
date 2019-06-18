/* History: Aug 13 2018 Initial coding
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define DEBUG 0

#define CHECK_MEM(obj) if (obj == NULL) {error("Memory");}
#define TWOPI 6.2831853071795862
#define LOG2 0.69314718055994529
#define LOG0ARG 1.0e-300
#define LOG0 -690.77552789821368
#define MAX_EP5 1.0e300

/* For optimizer */
#define stepredn	0.2
#define acctol	0.0001
#define reltest	10.0
#define expTo0       -9.0e100

/*
void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
*/


/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }

  return(ret);

} /* END: dVec_alloc */

/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate a double array (3-dim) */
/*
static double *** dArray_alloc(n1, n2, n3, initFlag, initVal)
int n1, n2, n3, initFlag;
double initVal;
{
  double ***mat;
  int i, j;

  mat = (double ***) malloc(n1*sizeof(double **));
  CHECK_MEM(mat);
  for (i=0; i<n1; i++) {
    mat[i] = (double **) malloc(n2*sizeof(double *));
    CHECK_MEM(mat[i]);
    for (j=0; j<n2; j++) mat[i][j] = dVec_alloc(n3, initFlag, initVal);
  }
     
  return(mat);

}  */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to free an aray (3d) */
/*
static void array_free(x, n1, n2)
void ***x;
int n1, n2;
{
  int i, j;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2; j++) free(x[i][j]);
    free(x[i]);
  }
  free(x);

}  */


/* Function to fill in a matrix from a vector (by column) */
static void fillMat(vec, nr, nc, out)
double *vec, **out;
int nr, nc;
{
  int i, j, col=0, ii;

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++;
    }
    col++;
  }

} /* END: fillMat */

 


static void get_pzzks(n, J, hatp0, akh_all, bkh_all, tmp0, tmp1, tmp2, tmp3, 
             t1, t2, lr, logor2, mu0_lrs, notMissing, ncpVecJ, vec, tvec, ret)
int n, J, *notMissing;
double t1, t2, tmp0, tmp1, tmp2, tmp3, **hatp0, **akh_all, **bkh_all;
double *logor2, *lr, *mu0_lrs, *ncpVecJ, *vec, *tvec, **ret;
{
  int k, i, j;
  double *vd1, xchisq, lrk, d1, d2, *bvec, *MAT, *avec, *pd;

  vd1 = dVec_alloc(J, 0, 0.0); 
  MAT = dVec_alloc(J, 0, 0.0); 

  for (k=1; k<n-1; k++) {
    /* Get tmat[k, ] */
    if (notMissing[k]) {
      xchisq = t2*logor2[k];
      for (i=0; i<J; i++) {
        vd1[i] = t2*dnchisq(xchisq, 1, ncpVecJ[i], 0);
      }
    } else {
      for (i=0; i<J; i++) vd1[i] = 1.0;
    }

    lrk  = lr[k];
    bvec = bkh_all[k];
    for (i=0; i<J; i++) {
      d1     = lrk-mu0_lrs[i];
      d2     = d1*d1;
      MAT[i] = bvec[i]*(tmp1/(tmp2*pow(tmp0+d2/t1, tmp3))*vd1[i]);
    }

    /* pzzks[,,k-1] <- tvec[k-1]*(hatp0*(akh_all[k-1,]%*%mat)) */

    d1  = tvec[k-1];
    pd = akh_all[k-1]; 
    for (i=0; i<J; i++) {
      d2 = pd[i];
      for (j=0; j<J; j++) {
        ret[i][j] += d1*hatp0[i][j]*d2*MAT[j];
      }
    }
  }


  avec = akh_all[n-2];
  for (i=0; i<J; i++) {
    d1 = avec[i];
    for (j=0; j<J; j++) {
      ret[i][j] += hatp0[i][j]*d1*vec[j];
    }
  }

  free(vd1);
  free(MAT);

} /* END: get_pzzks */

static double sum_pzzks(pzzks, n11, n12, n21, n22)
double **pzzks;
int n11, n12, n21, n22;
{
  int i, j;
  double sum=0.0;

  for (i=n11; i<=n12; i++) {
    for (j=n21; j<=n22; j++) sum += pzzks[i][j];
    
  }

  return(sum);

} /* END: sum_pzzks */

/* Function to compute 
  hatp0_m  <- c(1-hatp0_m01, hatp0_m01, hatp0_m10, 1-hatp0_m10)
*/
static void get_hatp0_m(pzzks, n, J, mainJ, ret)
double **pzzks, *ret;
int n, J, mainJ;
{
  /* 
    hatp0_m01 <- sum(pzzks[1:mainJ,(mainJ+1):J,])/sum(pzzks[1:mainJ,,])
    hatp0_m10 <- sum(pzzks[(mainJ+1):J, 1:mainJ,])/sum(pzzks[(mainJ+1):J,,])
    hatp0_m   <- c(1-hatp0_m01, hatp0_m01, hatp0_m10, 1-hatp0_m10)
  */
 
  double d1, d2, d3, d4, hatp0_m01, hatp0_m10;

  d1 = sum_pzzks(pzzks, 0,     mainJ-1, mainJ, J-1);
  d2 = sum_pzzks(pzzks, 0,     mainJ-1, 0,     J-1);
  d3 = sum_pzzks(pzzks, mainJ, J-1,     0,     mainJ-1);
  d4 = sum_pzzks(pzzks, mainJ, J-1,     0,     J-1);

  hatp0_m01 = d1/d2;
  hatp0_m10 = d3/d4;

  ret[0] = 1-hatp0_m01;
  ret[1] = hatp0_m01;
  ret[2] = hatp0_m10;
  ret[3] = 1-hatp0_m10;

  return;

} /* END: get_hatp0_m */

static void get_hatp0_M(pzzks, n, J, mainJ, tauJ, ret)
double **pzzks, *ret;
int n, J, mainJ, tauJ;
{
  int i1, g, l, i2, i3, i41, i42, i51, i52;
  double num, den, d1, d2, d3, d4, d5, d6, *pret;

  /* ret here is a vector but will be converted to an mainJ by tauJ
     matrix in the R code. In the R code, it is filled in by column,
     same below */

  pret = ret;
  i1   = 1 + mainJ;
  for (g=1; g<=mainJ; g++) {
    i2  = (g-1)*tauJ;
    i3  = g*tauJ;
    i41 = i1 + i2;
    i42 = i1 + i3 - 1;
    for (l=1; l<=mainJ; l++) {
      i51 = i1 + (l-1)*tauJ;
      i52 = i1 + l*tauJ - 1;

      /* Offset array indices by 1 */
      d1  = sum_pzzks(pzzks, l-1,   l-1,   g-1,   g-1);
      d2  = sum_pzzks(pzzks, l-1,   l-1,   i41-1, i42-1);
      d3  = sum_pzzks(pzzks, i51-1, i52-1, g-1,   g-1);
      d4  = sum_pzzks(pzzks, i51-1, i52-1, i41-1, i42-1);
      d5  = sum_pzzks(pzzks, l-1,   l-1,   0,     J-1);
      d6  = sum_pzzks(pzzks, i51-1, i52-1, 0,     J-1);
      num = d1 + d2 + d3 + d4;
      den = d5 + d6; 
      *pret++ = num/den;
    }
  }

  return;

} /* END: get_hatp0_M */

static double sum_pzks1(pzks1, n11, n12, n21, n22)
double **pzks1;
int n11, n12, n21, n22;
{
  int i, j;
  double sum=0.0;

  for (i=n11; i<=n12; i++) {
    for (j=n21; j<=n22; j++) sum += pzks1[i][j];
  }

  return(sum);

} /* sum_pzks1 */

static void get_r0s(pzks1, n, J, mainJ, tauJ, hatp0_mgtype, ret)
double **pzks1, **hatp0_mgtype, *ret;
int n, J, mainJ, tauJ;
{
  double d0, d1, d2, *r0, u01, u02;
  int i, j, i21, i22, ii;

  r0  = dVec_alloc(mainJ, 0, 0.0);
  d0  = sum_pzks1(pzks1, 0, 0, 0, J-1);

  /* Compute r0 */
  for (j=1; j<=mainJ; j++) {
    d1      = pzks1[0][j-1];
    i21     = mainJ+1+tauJ*(j-1);
    i22     = mainJ+tauJ*j; 
    d2      = sum_pzks1(pzks1, 0, 0, i21-1, i22-1);
    r0[j-1] = (d1 + d2)/d0; 
  }

  /* u0 <- c(sum(pzks1[1,1:mainJ]), sum(pzks1[1,(mainJ+1):J]))/sum(pzks1[1,]) */
  u01 = sum_pzks1(pzks1, 0, 0, 0, mainJ-1)/d0;
  u02 = sum_pzks1(pzks1, 0, 0, mainJ, J-1)/d0;

  /* 
   r0s  <- matrix(rep(r0*u0[2], each=tauJ), nrow=mainJ, ncol=tauJ, byrow=TRUE)*hatp0_mgtype
   r0s  <- c(r0*u0[1], as.vector(t(r0s)))
  */

  for (j=0; j<mainJ; j++) ret[j] = r0[j]*u01;
  ii = mainJ;
  for (i=0; i<mainJ; i++) {
    d1 = r0[i]*u02;
    for (j=0; j<tauJ; j++) {
      ret[ii] = d1*hatp0_mgtype[i][j];
      ii++;
    }
  }


  free(r0);

  return;

} /* END: get_r0s */

static void get_hatp0_mgtype(pzzks, n, J, mainJ, tauJ, r0s, ret)
double **pzzks, *r0s, *ret;
int n, J, mainJ, tauJ;
{
  double d2, d3, den, *pret;
  int i, j, i21, i22, ii;

  pret = ret;
  for (j=1; j<=mainJ; j++) {
    i21 = 1+mainJ+tauJ*(j-1);
    i22 = mainJ+tauJ*j;
    d3  = sum_pzzks(pzzks, 0, J-1, i21-1, i22-1);
    d2  = 0.0;
    for (i=i21-1; i<i22; i++) d2 += r0s[i];
    den = d2 + d3;
    ii  = i21-1;
    for (i=0; i<tauJ; i++) {
      d3      = sum_pzzks(pzzks, 0, J-1, ii, ii);
      *pret++ = (r0s[ii] + d3)/den;
      ii++;
    }
  }

  return;

} /* END: get_hatp0_mgtype */

void locFunc_hatp0(pn, pJ, pmainJ, ptauJ, pt1, pt2, tvec, ncpVecJ,
                   ptmp0, ptmp1, ptmp2, ptmp3, logor2, lr, notMissing,
                   akh_allVec, bkh_allVec, mu0_lrs, hatp0Vec, vec, 
                   hatp0_mgtypeVec, pzks1Vec,
                   ret_code, ret_hatp0_M, ret_hatp0_m, ret_r0s, ret_hatp0_mgtype)
int *pn, *pJ, *pmainJ, *ptauJ, *notMissing, *ret_code;
double *pt1, *pt2, *tvec, *ncpVecJ, *ptmp0, *ptmp1, *ptmp2, *ptmp3;
double *logor2, *lr, *akh_allVec, *bkh_allVec, *mu0_lrs, *hatp0Vec, *vec;
double *hatp0_mgtypeVec, *pzks1Vec;
double *ret_hatp0_M, *ret_hatp0_m, *ret_r0s, *ret_hatp0_mgtype;
{
  int n, J, mainJ, tauJ;
  double t1, t2, tmp0, tmp1, tmp2, tmp3;
  double **akh_all, **bkh_all, **hatp0, **hatp0_mgtype, **pzks1;
  double **pzzks;

  *ret_code = -1;
  n         = *pn;
  J         = *pJ;
  mainJ     = *pmainJ;
  tauJ      = *ptauJ;
  t1        = *pt1;
  t2        = *pt2;
  tmp0      = *ptmp0;
  tmp1      = *ptmp1;
  tmp2      = *ptmp2;
  tmp3      = *ptmp3;

  akh_all      = dMat_alloc(n, J, 0, 0.0);
  bkh_all      = dMat_alloc(n, J, 0, 0.0);
  hatp0        = dMat_alloc(J, J, 0, 0.0);
  hatp0_mgtype = dMat_alloc(mainJ, tauJ, 0, 0.0);
  pzks1        = dMat_alloc(n, J, 0, 0.0);
  pzzks        = dMat_alloc(J, J, 1, 0.0);

  fillMat(akh_allVec, n, J, akh_all);
  fillMat(bkh_allVec, n, J, bkh_all);
  fillMat(hatp0Vec, J, J, hatp0);
  fillMat(hatp0_mgtypeVec, mainJ, tauJ, hatp0_mgtype);
  fillMat(pzks1Vec, n, J, pzks1);


  get_pzzks(n, J, hatp0, akh_all, bkh_all, tmp0, tmp1, tmp2, tmp3, 
    t1, t2, lr, logor2, mu0_lrs, notMissing, ncpVecJ, vec, tvec, pzzks);
  get_hatp0_m(pzzks, n, J, mainJ, ret_hatp0_m);
  get_hatp0_M(pzzks, n, J, mainJ, tauJ, ret_hatp0_M);
  get_r0s(pzks1, n, J, mainJ, tauJ, hatp0_mgtype, ret_r0s);
  get_hatp0_mgtype(pzzks, n, J, mainJ, tauJ, ret_r0s, ret_hatp0_mgtype);


  
  matrix_free((void **) akh_all, n); 
  matrix_free((void **) bkh_all, n);
  matrix_free((void **) hatp0, J);
  matrix_free((void **) hatp0_mgtype, mainJ);
  matrix_free((void **) pzks1, n);
  matrix_free((void **) pzzks, J, J);


  *ret_code = 0;
  return;

} /* END: locFunc_hatp0 */

static void rowVecMat(vec, mat, nr, nc, ret)
double *vec, **mat, *ret;
int nr, nc;
{
  int i, j;
  double sum;

  for (i=0; i<nc; i++) {
    sum = 0.0;
    for (j=0; j<nr; j++) {
      sum += vec[j]*mat[j][i];
    }
    ret[i] = sum;
  }

  return;

} /* END: rowVecMat */

static void matColVec(vec, mat, nr, nc, ret)
double *vec, **mat, *ret;
int nr, nc;
{
  int i, j;
  double sum;

  for (i=0; i<nr; i++) {
    sum = 0.0;
    for (j=0; j<nc; j++) {
      sum += vec[j]*mat[i][j];
    }
    ret[i] = sum;
  }

  return;

} /* END: matColVec */


void akh_ca(pn, pJ, pmainJ, ptauJ, 
                   ppt5v0, ppt5v0POWpt5v0, pv0Plus1Over2, pgamma_v0Over2,
                   ptmp1, ptmp2, ptmp3, tmp4, logor2, lr, notMissing,
                   mu0_lrs, a1h, hatp0Vec, pc1, 
                   ret_code, ret_akh_all, ret_ca_all)
int *pn, *pJ, *pmainJ, *ptauJ, *notMissing, *ret_code;
double *ppt5v0, *ppt5v0POWpt5v0, *pv0Plus1Over2, *pgamma_v0Over2;
double *ptmp1, *ptmp2, *ptmp3, *tmp4;
double *logor2, *lr, *a1h, *mu0_lrs, *hatp0Vec, *pc1;
double *ret_akh_all, *ret_ca_all;
{
  int n, J, k, ii, j, notMiss;
  double tmp1, tmp2, tmp3, c1, tmp5, cka, lrk, logor2k, sum_akd;
  double **hatp0, pt5v0, v0Plus1Over2, gamma_v0Over2;
  double *a1h_hatp0, *akd, d1, d2;

  *ret_code     = -1;
  n             = *pn;
  J             = *pJ;
  tmp1          = *ptmp1;
  tmp2          = *ptmp2;
  tmp3          = *ptmp3;
  c1            = *pc1;
  pt5v0         = *ppt5v0;
  v0Plus1Over2  = *pv0Plus1Over2;
  gamma_v0Over2 = *pgamma_v0Over2;

  /* akh_all.new[1,] <- a1h, ca_all.new[1]   <- c1 */
  ret_ca_all[0] = c1;
  for (k=0; k<J; k++) ret_akh_all[k] = a1h[k];

  hatp0     = dMat_alloc(J, J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);
  a1h_hatp0 = dVec_alloc(J, 0, 0.0);
  akd       = dVec_alloc(J, 0, 0.0);

  /* Initialize akd with a1h */
  for (j=0; j<J; j++) akd[j] = a1h[j];

  ii = J;
  for (k=1; k<n; k++) {
    notMiss = notMissing[k];
    lrk     = lr[k];
    logor2k = logor2[k];
    sum_akd = 0.0;

    /* Compute a1h*hatp0 */
    rowVecMat(akd, hatp0, J, J, a1h_hatp0);

    for (j=0; j<J; j++) {
      d1   = lrk-mu0_lrs[j];
      d1   = d1*d1;
      d2   = pow(pt5v0+(d1/tmp2), v0Plus1Over2);
      tmp5 = tmp1/(gamma_v0Over2*d2);
      
      if (!notMiss) {
        d1 = a1h_hatp0[j]*tmp5;
      } else { 
        d1 = a1h_hatp0[j]*(tmp5*(tmp3*dnchisq(logor2k*tmp3, 1, tmp4[j], 0)));
      }
      sum_akd += d1;
      akd[j]  = d1;
    }
    cka             = 1.0/sum_akd;
    ret_ca_all[k]   = cka;
    for  (j=0; j<J; j++) {
      d1              = cka*akd[j];
      akd[j]          = d1;
      ret_akh_all[ii] = d1;
      ii++;
    }
  }

  matrix_free((void **) hatp0, J);
  free(a1h_hatp0);
  free(akd);

  *ret_code = 0;
  return;

} /* END: akh_ca */

void C_hes_logL1(pn, pJ, parm, mu0_lrs, mu0_lors,
               ppi, logor2, lr, notMissing, r0s, hatp0Vec,  
               ret_code, ret_loglike)
int *pn, *pJ, *notMissing, *ret_code;
double *parm, *mu0_lrs, *mu0_lors, *ppi, *logor2, *lr, *r0s, *hatp0Vec; 
double *ret_loglike;
{
  int n, J, k, j, notMiss;
  double lrk, logor2k, sum_akd, dpi, lsigma20_lr, lsigma20_lor, v0, sigma20_lr;
  double **hatp0, p5v0, p5v0p5v0, v01, gv01, gv0, tmp1, tmp2, tmp3, cka, *a1d;
  double *a1h_hatp0, sum, akd, d1, d2, d3, *ncpVec, logL_Y, c1, sum_log_ca_all, dtmp;

#if DEBUG
Rprintf("BEGIN: C_hes_logL1\n");
#endif

  *ret_code     = -1;
  n             = *pn;
  J             = *pJ;
  dpi           = *ppi;
  lsigma20_lr   = parm[3-1];
  lsigma20_lor  = parm[4-1];
  v0            = parm[5-1];
  sigma20_lr    = exp(lsigma20_lr);
  p5v0          = 0.5*v0;
  p5v0p5v0      = pow(p5v0, p5v0);
  v01           = (v0 + 1.0)/2.0;
  gv01          = gammafn(v01);
  gv0           = gammafn(v0/2.0);
  tmp1          = p5v0p5v0*gv01*pow(2.0*dpi*sigma20_lr, -0.5);
  tmp2          = 2.0*sigma20_lr;
  tmp3          = exp(-lsigma20_lor);

  ncpVec    = dVec_alloc(J, 0, 0.0);
  a1d       = dVec_alloc(J, 0, 0.0);
  a1h_hatp0 = dVec_alloc(J, 0, 0.0);
  hatp0     = dMat_alloc(J, J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);

#if DEBUG
Rprintf("Compute ncp\n");
#endif

  for (k=0; k<J; k++) ncpVec[k] = mu0_lors[k]*mu0_lors[k]*tmp3*r0s[k];

#if DEBUG
Rprintf("Compute a1d\n");
#endif

  /* Compute a1d */
  notMiss = notMissing[0];
  lrk     = lr[0];
  logor2k = logor2[0];
  sum     = 0.0;
  for (k=0; k<J; k++) {
    d1 = lrk - mu0_lrs[k];
    d1 = d1*d1;
    d2 = gv0*pow(p5v0 + d1/tmp2, v01);
    d3 = (tmp1/d2)*r0s[k];  
    if (notMiss) d3 = d3*tmp3*dnchisq(logor2k*tmp3, 1.0, ncpVec[k], 0);

    a1d[k] = d3;
    sum += d3;
  }
/*
print_dVec(a1d, J, "a1d");
*/


#if DEBUG
Rprintf("Re-define ncp\n");
#endif

  c1             = 1.0/sum;
  sum_log_ca_all = log(c1);
  for (k=0; k<J; k++) {
    a1d[k]    = c1*a1d[k]; 
    ncpVec[k] = mu0_lors[k]*mu0_lors[k]*tmp3;
  }

#if DEBUG
dtmp = ncpVec[0];
for (k=1; k<J; k++) {
  d1 = ncpVec[k];
  if (d1 > dtmp) dtmp = d1;
}
Rprintf("Max ncp = %g\n", dtmp);
Rprintf("Main loop\n");
#endif


  for (k=1; k<n; k++) {
    notMiss = notMissing[k];
    logor2k = logor2[k];
    lrk     = lr[k];
    dtmp    = logor2k*tmp3;
    
    rowVecMat(a1d, hatp0, J, J, a1h_hatp0);
    sum_akd = 0.0;
    for (j=0; j<J; j++) {
      d1  = lrk - mu0_lrs[j];
      d1  = d1*d1;
      d2  = gv0*pow(p5v0 + d1/tmp2, v01);
      akd = a1h_hatp0[j]*(tmp1/d2);  
      if (notMiss) akd = akd*tmp3*dnchisq(dtmp, 1.0, ncpVec[j], 0);
      sum_akd += akd;
      a1d[j] = akd;
    }
    cka = 1.0/sum_akd;

    sum_log_ca_all += log(cka);
    sum_akd = 0.0;
    for (j=0; j<J; j++) {
      d1       = a1d[j]*cka;
      a1d[j]   = d1;
      sum_akd += d1;
    }
  }

#if DEBUG
Rprintf("end main loop\n");
#endif


  logL_Y = -sum_log_ca_all + log(sum_akd);
  *ret_loglike = logL_Y;
  
  free(ncpVec);
  free(a1d);
  free(a1h_hatp0);
  matrix_free((void **) hatp0, J);

  *ret_code = 0;

#if DEBUG
Rprintf("END: C_hes_logL1\n");
#endif


  return;

} /* END: C_hes_logL1 */

void bkh_cb(pn, pJ, ptmp0, ptmp1, ptmp2, ptmp3, ptmp4, ptmp5, tmp6,
    logor2, lr, notMissing, mu0_lrs, hatp0Vec,  
                   ret_code, ret_bkh_all, ret_cb_all)
int *pn, *pJ, *notMissing, *ret_code;
double *ptmp0, *ptmp1, *ptmp2, *ptmp3, *ptmp4, *ptmp5, *tmp6;
double *logor2, *lr, *mu0_lrs, *hatp0Vec;
double *ret_bkh_all, *ret_cb_all;
{
  int n, J, k, ii, kp1, j;
  double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, lrk1, logor2k1;
  double **hatp0, d1, d2, *bnh, *vec, *bkd, ckb;
  
  *ret_code = -1;
  n         = *pn;
  J         = *pJ;
  tmp0      = *ptmp0;
  tmp1      = *ptmp1;
  tmp2      = *ptmp2;
  tmp3      = *ptmp3;
  tmp4      = *ptmp4;
  tmp5      = *ptmp5;

  hatp0     = dMat_alloc(J, J, 0, 0.0);
  bnh       = dVec_alloc(J, 0, 0.0);
  vec       = dVec_alloc(J, 0, 0.0);
  bkd       = dVec_alloc(J, 0, 0.0);
  fillMat(hatp0Vec, J, J, hatp0);

  /* We are working backwards, so get results for final index, row */
  d1 = 1.0/((double) J);
  ret_cb_all[n-1] = d1;
  ii = (n-1)*J;
  for (k=0; k<J; k++) {
    ret_bkh_all[ii] = d1;
    ii++;
    bnh[k] = d1;
  }

  for (k=n-2; k>-1; k--) {
    kp1      = k + 1;
    logor2k1 = logor2[kp1];
    lrk1     = lr[kp1];

    if (!notMissing[kp1]) {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*bnh[j];
      }
    } else {
      for (j=0; j<J; j++) {
        d1     = lrk1 - mu0_lrs[j];
        d2     = pow(tmp0 + d1*d1/tmp3, tmp4);
        vec[j] = tmp1/(tmp2*d2)*(tmp5*dnchisq(logor2k1*tmp5, 1.0, tmp6[j], 0))*bnh[j];
      }
    }
    matColVec(vec, hatp0, J, J, bkd);
    d1  = 0.0;
    for (j=0; j<J; j++) d1 += bkd[j];
    ckb = 1.0/d1;

    ret_cb_all[k] = ckb;
    ii            = k*J;   
    for (j=0; j<J; j++) {
      d1 = ckb*bkd[j];
      ret_bkh_all[ii] = d1;
      ii++;
      bnh[j] = d1;
    }
  } 

  free(bkd);
  free(bnh);
  free(vec);
  matrix_free((void **) hatp0, J);

  *ret_code = 0;
  return;
  
} /* END: bkh_cb */



static double ** Lmatrix(int n)
{
    int   i;
    double **m;

    m = (double **) malloc(n*sizeof(double *));
    for (i = 0; i < n; i++)
	m[i] = (double *) malloc((i + 1)*sizeof(double));
    return m;
}





static double logBase2(d)
double d;
{
  double ret;
  ret = log(d)/LOG2;
  return(ret);

} /* END: logBase2 */

static void mu0_lr_func(nj, mainJ, lpsi0, lpr0, lsz0, ctzs, ctms, ret)
int nj, mainJ;
double lpsi0, lpr0, lsz0, *ctzs, *ctms, *ret;
{
  int i;
  double pr0, sz0, tmp1;
  
  pr0 = 1.0/(1.0 + exp(lpr0));
  sz0 = 1.0/(1.0 + exp(lsz0));
  for (i=0; i<mainJ; i++) ret[i] = logBase2(2.0+pr0*(ctzs[i]-2.0)) - lpsi0;
  for (i=mainJ; i<nj; i++) {
    tmp1   = ctms[i - mainJ];
    ret[i] = logBase2(2.0+(tmp1-2.0)*pr0+(ctzs[i]-tmp1)*sz0*pr0) - lpsi0;
  } 
  
  return;

} /* END: mu0_lr_func */    

static void mu0_lor_func(nj, mainJ, lpr0, lsz0, mz, pz, mz_sub, pz_sub, ret)
int nj, mainJ;
double lpr0, lsz0, *mz, *pz, *mz_sub, *pz_sub, *ret;
{
  int i, j;
  double pr0, sz0, tmp1, tmp2, arg1, arg2, oneMinusPr0, szpr;
  
  pr0         = 1.0/(1.0 + exp(lpr0));
  sz0         = 1.0/(1.0 + exp(lsz0));
  oneMinusPr0 = 1.0 - pr0;
  szpr        = sz0*pr0;

  for (i=0; i<mainJ; i++) {
    arg1 = mz[i]*pr0 + oneMinusPr0;
    if (arg1 < LOG0ARG) {
      tmp1 = LOG0; 
    } else {
      tmp1 = log(arg1);
    }
    arg2 = pz[i]*pr0 + oneMinusPr0;
    if (arg2 < LOG0ARG) {
      tmp2 = LOG0; 
    } else {
      tmp2 = log(arg2);
    }
    ret[i] = tmp1 - tmp2;
  }

  for (i=mainJ; i<nj; i++) {
    j    = i - mainJ;
    tmp1 = mz_sub[j];
    arg1 = 1.0 + (mz[i] - tmp1)*szpr + (tmp1 - 1.0)*pr0;   
    if (arg1 < LOG0ARG) {
      tmp1 = LOG0; 
    } else {
      tmp1 = log(arg1);
    }
    tmp2 = pz_sub[j];
    arg2 = 1.0 + (pz[i] - tmp2)*szpr + (tmp2 - 1.0)*pr0;
    if (arg2 < LOG0ARG) {
      tmp2 = LOG0; 
    } else {
      tmp2 = log(arg2);
    }
    ret[i] = tmp1 - tmp2;
  } 
  
  return;

} /* END: mu0_lor_func */   

static double negLoglike(parm, nparm, p5L, p6L, p6U, n, J, mainJ, lr, logor2,
    notMissing, pzks1, ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors)
int nparm, n, J, *notMissing, mainJ;
double *parm, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub;
double *mu0_lrs, *mu0_lors;
{
  int i, nomiss, j;
  double p1, p2, p3, p4, p5, p6, p62, sigma20_lr, ep5, d1, d2, d3, mat, ret=0.0, p6p1, lri;
  double *vec, pzks1ij, logor2i, darg;

  p1         = parm[0];
  p2         = parm[1];
  p3         = parm[2];
  p4         = parm[3];
  p5         = parm[4];
  p6         = parm[5];
  sigma20_lr = exp(p4);
  
  /* Reparameterize */
  p62  = p6*p6;
  p6   = p6L + (p6U - p6L)*p62/(1.0 + p62);

  p5   = p5L + p5*p5;
  ep5  = exp(-p5);
  if (ep5 > MAX_EP5) error("ERROR: set bounds on parm[5]");
  p6p1 = p6 + 1.0;

  mu0_lr_func(J, mainJ, p1, p2, p3, ctzs, ctms, mu0_lrs);
  mu0_lor_func(J, mainJ, p2, p3, mz, pz, mz_sub, pz_sub, mu0_lors);


  /* Define the ncp and store in mu0_lors */
  for (i=0; i<J; i++) {
    d1          = mu0_lors[i];
    mu0_lors[i] = d1*d1*ep5;
  }

  d1 = -0.5*log(TWOPI*sigma20_lr) + log(gammafn(p6p1/2.0)) - log(gammafn(p6/2.0)) + 0.5*p6*log(p6/2.0);
  d2 = 0.5*p6p1;
  for (i=0; i<n; i++) {
    lri     = lr[i];
    logor2i = logor2[i];
    vec     = pzks1[i];
    nomiss  = notMissing[i];
    darg    = logor2i*ep5;
    for (j=0; j<J; j++) {
      d3      = lri - mu0_lrs[j];
      pzks1ij = vec[j];
      mat     = (d1 - d2*log(0.5*(p6 + d3*d3/sigma20_lr)))*pzks1ij;
      if (nomiss) mat = mat - pzks1ij*(p5 - dnchisq(darg, 1.0, mu0_lors[j], 1));
      ret += mat;
    }
  }
  ret = -ret;

  return(ret);

} /* END: negLoglike */

static void gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, ret)
double *eta, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub, *mu0_lrs, *mu0_lors;
int nparm, nobs, J, mainJ, *notMissing;
double *ret;
{
  int i;
  double fxplush, fxminush, h, save, h2, *ptr, *ptrret;

  /* default step size */
  h = 1e-3;
  h2 = 2.0*h;

  for (i=0, ptr=eta, ptrret=ret; i<nparm; i++, ptr++, ptrret++) {
    save = *ptr;  

    *ptr = save + h;  
    /*fxplush = negloglike(eta, op);*/
    fxplush = negLoglike(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

    *ptr = save - h;  
    /*fxminush = negloglike(eta, op);*/
    fxminush =  negLoglike(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
    *ptr = save;

    *ptrret = (fxplush - fxminush)/h2;
  }

} /* END: gradient */


static void myvmmin(int n0, double *b, double *Fmin, 
      int maxit, int trace, double abstol, double reltol, int nREPORT, 
      int *fncount, int *grcount, int *fail,
      double p5L, double p6L, double p6U, int nobs, int J, int mainJ, double *lr, double *logor2, int *notMissing, 
      double **pzks1, double *ctzs, double *ctms, double *mz, double *pz, 
      double *mz_sub, double *pz_sub, double *mu0_lrs, double *mu0_lors)
{
    int accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l, nparm=n0;

    if (maxit <= 0) {
	*fail = 0;
	/* *Fmin = negloglike(b, op);*/
       *Fmin =  negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
	*fncount = *grcount = 0;
	return;
    }

    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    /*for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;*/
    for (i = 0; i < n0; i++) l[n++] = i;
    g = (double *) R_alloc(n0, sizeof(double));
    t = (double *) R_alloc(n, sizeof(double));
    X = (double *) R_alloc(n, sizeof(double));
    c = (double *) R_alloc(n, sizeof(double));
    B = Lmatrix(n);
    /*f = negloglike(b, op);*/
    f = negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);
 
    if (!isfinite(f))
	Rprintf("initial value in 'vmmin' is not finite");
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    /*gradient(b, op, g);*/
    gradient(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = 0;
	    do {

		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
		    /*f = negloglike(b, op);*/
                  f = negLoglike(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
                              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

		    funcount++;
		    accpoint = isfinite(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (!accpoint) {
			steplength *= stepredn;
		    }
		}

	    } while (!(count == n || accpoint));
	    enough = (f > abstol) && 
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);

	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }

	    if (count < n) {/* making progress */
		*Fmin = f;
		/*gradient(b, op, g);*/
              gradient(b, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, g);

		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
    /*
    if (op->debug) {
      Rprintf("fncount = %d   grcount = %d \n", funcount, gradcount);
    }
    */
}


/* Main function for optimizing loglike */
void myoptimC(parm, pnparm, pp5L, pp6L, pp6U, pnobs, pJ, pmainJ, lr, logor2, notMissing,
             pzks1Vec, ctzs, ctms, mz, pz, mz_sub, pz_sub, preltol,
             ret_code)
double *parm, *pp5L, *pp6L, *pp6U, *lr, *logor2, *pzks1Vec, *ctzs, *ctms;
double *mz, *pz, *mz_sub, *pz_sub, *preltol;
int *pnparm, *pJ,*pnobs, *pmainJ, *notMissing, *ret_code; 
{
  int nparm, nobs, J, mainJ;
  int nREPORT, fail, maxit, trace, fncount, grcount;
  double p5L, p6L, p6U, abstol, reltol, Fmin, **pzks1, *mu0_lrs, *mu0_lors;

  *ret_code = -1;
  nparm     = *pnparm;
  nobs      = *pnobs;
  p5L       = *pp5L;
  p6L       = *pp6L;
  p6U       = *pp6U;
  J         = *pJ;
  mainJ     = *pmainJ;
  reltol    = *preltol;

  nREPORT   = 1;
  abstol    = -1.0e100;
  fail      = 1;
  fncount   = 0;
  grcount   = 0;
  maxit     = 100;
  trace     = 0;

  mu0_lrs  = dVec_alloc(J, 0, 0.0);
  mu0_lors = dVec_alloc(J, 0, 0.0);
  pzks1    = dMat_alloc(nobs, J, 0, 0.0);
  fillMat(pzks1Vec, nobs, J, pzks1);

  myvmmin(nparm, parm, &Fmin, maxit, trace, abstol, reltol, nREPORT, 
      &fncount, &grcount, &fail,
      p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, 
      pzks1, ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors);

  free(mu0_lrs);
  free(mu0_lors);
  matrix_free((void **) pzks1, nobs);

  *ret_code = fail;

  return;

} /* END: myoptimC */

/*
mu0_lr_func(nj, mainJ, lpsi0, lpr0, lsz0, ctzs, ctms, ret)
mu0_lor_func(nj, mainJ, lpr0, lsz0, mz, pz, mz_sub, pz_sub, ret)

hes_logL1(pn, pJ, parm, mu0_lrs, mu0_lors,
               ppi, logor2, lr, notMissing, r0s, hatp0Vec,  
               ret_code, ret_loglike)


static void optim_hessian(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, ret)
double *eta, p5L, p6L, p6U, *lr, *logor2, **pzks1, *ctzs, *ctms, *mz, *pz, *mz_sub, *pz_sub, *mu0_lrs, *mu0_lors;
int nparm, nobs, J, mainJ, *notMissing;
double *ret;
{
  int i, nparms, j, ni;
  double *gxplush, *gxminush, h, save, h2, *ptr;

  gxplush  = dVec_alloc(nparms, 0, 0.0);
  gxminush = dVec_alloc(nparms, 0, 0.0);

  h = 1e-3;
  h2 = 2.0*h;
  for (i=0, ptr=eta; i<nparms; i++, ptr++) {
    save = *ptr;  

    *ptr = save + h;  
    gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, gxplush);

    *ptr = save - h;  
    gradient(eta, nparm, p5L, p6L, p6U, nobs, J, mainJ, lr, logor2, notMissing, pzks1, 
              ctzs, ctms, mz, pz, mz_sub, pz_sub, mu0_lrs, mu0_lors, gxminush);

    *ptr = save;

    ni = i*nparms;
    for (j=i; j<nparms; j++) ret[ni + j] = (gxplush[j] - gxminush[j])/h2;
    for (j=i+1; j<nparms; j++) ret[j*nparms + i] = ret[ni + j];
  }

  free(gxplush);
  free(gxminush);

}  END: optim_hessian */


