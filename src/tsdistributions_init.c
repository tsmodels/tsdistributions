// Dummy file required so that useDynLib(tsdistributions, .registration=TRUE) doesn't fail on empty 'src'
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL
#include <R_ext/Visibility.h>
/* .C calls */
extern void c_dghst(double *, double *, double *, double *, double *, double *, int *, int*);
extern void c_rghst(int *, double *, double *, double *, double *, double *);
extern void c_dghyp(double *, double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_dgh(double *, double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_rgh(int *, double *, double *, double *, double *, double *, double *);
extern void c_dnig(double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_rnig(int *, double *, double *, double *, double *, double *);
extern void c_rstd(int *, double *, double *, double *, double *);
extern void c_dstd(double *, double *, double *, double *, double *, int *, int *);
extern void c_pstd(double *, double *, double *, double *, double *, int *);
extern void c_qstd(double *, double *, double *, double *, double *, int *);
extern void c_rsstd(int *, double *, double *, double *, double *, double *);
extern void c_dsstd(double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_psstd(double *, double *, double *, double *, double *, double *, int *);
extern void c_qsstd(double *, double *, double *, double *, double *, double *, int *);
extern void c_djsu(double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_qjsu(double *, double *, double *, double *, double *, double *, int *);
extern void c_rjsu(int *, double *, double *, double *, double *, double *);
extern void c_pjsu(double *, double *, double *, double *, double *, double *, int *);
extern void c_rsnorm(int *, double *, double *, double *, double *);
extern void c_dsnorm(double *, double *, double *, double *, double *, int *, int *);
extern void c_psnorm(double *, double *, double *, double *, double *, int *);
extern void c_qsnorm(double *, double *, double *, double *, double *, int *);
extern void c_rged(int *, double *, double *, double *, double *);
extern void c_dged(double *, double *, double *, double *, double *, int *, int *);
extern void c_pged(double *, double *, double *, double *, double *, int *);
extern void c_qged(double *, double *, double *, double *, double *, int *);
extern void c_rsged(int *, double *, double *, double *, double *, double *);
extern void c_dsged(double *, double *, double *, double *, double *, double *, int *, int *);
extern void c_psged(double *, double *, double *, double *, double *, double *, int *);
extern void c_qsged(double *, double *, double *, double *, double *, double *, int *);
extern void c_dhyp(double *, double *, double *, double *, double *, double *, int *, int *);


static const R_CMethodDef CEntries[] = {
    {"c_dghst",          (DL_FUNC) &c_dghst,           8},
    {"c_rghst",          (DL_FUNC) &c_rghst,           6},
    {"c_dghyp",          (DL_FUNC) &c_dghyp,           9},
    {"c_dgh",            (DL_FUNC) &c_dgh,             9},
    {"c_dnig",           (DL_FUNC) &c_dnig,            8},
    {"c_rstd",           (DL_FUNC) &c_rstd,            5},
    {"c_dstd",           (DL_FUNC) &c_dstd,            7},
    {"c_pstd",           (DL_FUNC) &c_pstd,            6},
    {"c_qstd",           (DL_FUNC) &c_qstd,            6},
    {"c_rsstd",          (DL_FUNC) &c_rsstd,           6},
    {"c_dsstd",          (DL_FUNC) &c_dsstd,           8},
    {"c_psstd",          (DL_FUNC) &c_psstd,           7},
    {"c_qsstd",          (DL_FUNC) &c_qsstd,           7},
    {"c_djsu",           (DL_FUNC) &c_djsu,            8},
    {"c_qjsu",           (DL_FUNC) &c_qjsu,            7},
    {"c_rjsu",           (DL_FUNC) &c_rjsu,            6},
    {"c_pjsu",           (DL_FUNC) &c_pjsu,            7},
    {"c_rsnorm",         (DL_FUNC) &c_rsnorm,          5},
    {"c_dsnorm",         (DL_FUNC) &c_dsnorm,          7},
    {"c_psnorm",         (DL_FUNC) &c_psnorm,          6},
    {"c_qsnorm",         (DL_FUNC) &c_qsnorm,          6},
    {"c_rged",           (DL_FUNC) &c_rged,            5},
    {"c_dged",           (DL_FUNC) &c_dged,            7},
    {"c_pged",           (DL_FUNC) &c_pged,            6},
    {"c_qged",           (DL_FUNC) &c_qged,            6},
    {"c_rsged",          (DL_FUNC) &c_rsged,           6},
    {"c_dsged",          (DL_FUNC) &c_dsged,           8},
    {"c_psged",          (DL_FUNC) &c_psged,           7},
    {"c_qsged",          (DL_FUNC) &c_qsged,           7},
    {"c_dhyp",           (DL_FUNC) &c_dhyp,            8},
    {NULL, NULL, 0}
};

void attribute_visible R_init_tsdistributions(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
