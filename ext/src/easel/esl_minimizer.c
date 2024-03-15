/* Multidimensional optimization using conjugate gradient descent.
 * 
 * Can be used even without derivative information; falls back to
 * a numeric gradient if analytic gradient is unavailable.
 * 
 * Contents:
 *    1. esl_min_ConjugateGradientDescent(), our optimizer
 *    2. ESL_MIN_CFG, for optional configuration/customization
 *    3. ESL_MIN_DAT, for optional data collection
 *    4. Internal functions: numeric deriv, bracketing, 1D line min
 *    5. Unit tests
 *    6. Test driver
 *    7. Example
 */
#include <esl_config.h>

#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "esl_minimizer.h"



static void numeric_derivative(ESL_MIN_CFG *cfg, double *x, int n, 
			       double (*func)(double *, int, void*),
			       void *prm, double *dx, ESL_MIN_DAT *dat);
static int bracket(ESL_MIN_CFG *cfg, double *ori, double *d, int n, double firststep,
		   double (*func)(double *, int, void *), void *prm, double *wrk, 
		   double *ret_ax, double *ret_bx, double *ret_cx,
		   double *ret_fa, double *ret_fb, double *ret_fc, ESL_MIN_DAT *dat);
static void brent(ESL_MIN_CFG *cfg, double *ori, double *dir, int n,
		  double (*func)(double *, int, void *), void *prm, double a, double b, 
		  double *xvec, double *opt_x, double *opt_fx, ESL_MIN_DAT *dat);
 


/*****************************************************************
 * 1. esl_min_ConjugateGradientDescent(), our optimizer
 *****************************************************************/ 

/* Function:  esl_min_ConjugateGradientDescent()
 * Incept:    SRE, Wed Jun 22 08:49:42 2005 [St. Louis]
 *
 * Purpose:   n-dimensional minimization by conjugate gradient descent.
 *           
 *            An initial point is provided by <x>, a vector of <n>
 *            components. The caller also provides a function <*func()> that 
 *            compute the objective function f(x) when called as 
 *            <(*func)(x, n, prm)>, and a function <*dfunc()> that can
 *            compute the gradient <dx> at <x> when called as 
 *            <(*dfunc)(x, n, prm, dx)>, given an allocated vector <dx>
 *            to put the derivative in. Any additional data or fixed
 *            parameters that these functions require are passed by
 *            the void pointer <prm>.
 *            
 *            The first step of each iteration is to try to bracket
 *            the minimum along the current direction. The initial step
 *            size is controlled by <u[]>; the first step will not exceed 
 *            <u[i]> for any dimension <i>. (You can think of <u> as
 *            being the natural "units" to use along a graph axis, if
 *            you were plotting the objective function.)
 *
 *            The caller also provides an allocated workspace sufficient to
 *            hold four allocated n-vectors. (4 * sizeof(double) * n).
 *
 *            Iterations continue until the objective function has changed
 *            by less than a fraction <tol>. This should not be set to less than
 *            sqrt(<DBL_EPSILON>). 
 *
 *            Upon return, <x> is the minimum, and <opt_fx> (if
 *            provided) is f(x), the function value at <x>.
 *            
 * Args:      cfg      - optional custom config for tunable params; or NULL
 *            x        - an initial guess n-vector; RETURN: x at the minimum
 *            n        - dimensionality of all vectors
 *            *func()  - function for computing objective function f(x)
 *            *dfunc() - function for computing a gradient at x
 *            prm      - void ptr to any data/params func,dfunc need 
 *            opt_fx   - optRETURN: f(x) at the minimum; or NULL if unwanted
 *            opt_dat  - optRETURN: table of stats on the run; or NULL if unwanted
 *
 * Returns:   <eslOK> on success. <x> is updated to be the minimum. <*opt_fx>,
 *            if provided, contains f(x) at that minimum. <*opt_dat>, if
 *            an allocated <ESL_MIN_DAT> was provided, contains a table of 
 *            statistics on the run.
 *
 *            <eslENOHALT> if it fails to converge in max iterations,
 *            but the final <x>, <*opt_fx>, and <*opt_dat> are still provided;
 *            maybe they're good enough, caller can decide.
 *
 * Throws:    <eslERANGE> if the minimum is not finite, which may
 *            indicate a problem in the implementation or choice of <*func()>.
 *
 *            <eslEMEM> on allocation failure.
 *            On thrown exceptions, <*opt_fx> is eslINFINITY, and <x> is undefined.
 *
 * Xref:      STL9/101.
 */
int
esl_min_ConjugateGradientDescent(ESL_MIN_CFG *cfg, double *x, int n, 
       				 double (*func)(double *, int, void *),
				 void (*dfunc)(double *, int, void *, double *),
				 void *prm, double *opt_fx, ESL_MIN_DAT *dat)
{
  int     max_iterations = cfg ? cfg->max_iterations : eslMIN_MAXITER;
  double  cg_rtol        = cfg ? cfg->cg_rtol        : eslMIN_CG_RTOL;
  double  cg_atol        = cfg ? cfg->cg_atol        : eslMIN_CG_ATOL;
  double *u              = cfg ? cfg->u              : NULL;               // u[i] = 1.0 if not custom
  double *wrk            = NULL;
  double oldfx;
  double coeff;
  int    i, i1;
  double *dx, *cg, *w1, *w2;
  double fa,fb,fc;
  double ax,bx,cx;
  double fx;
  int    status;

  ESL_ALLOC(wrk, sizeof(double) * 4 * n);  // tmp workspace for 4 n-vectors
  dx = wrk;
  cg = wrk + n;
  w1 = wrk + 2*n;
  w2 = wrk + 3*n;

  oldfx = (*func)(x, n, prm);	/* init the objective function */
  if (dat) {  dat->fx[0] = oldfx; dat->nfunc[0] = 1; dat->niter = 0; }
  
  /* Bail out if the function is +/-inf or nan: this can happen if the caller
   * has screwed something up, or has chosen a bad start point.
   */
  if (! isfinite(oldfx)) ESL_XEXCEPTION(eslERANGE, "minimum not finite");

  if (dfunc)
    {
      (*dfunc)(x, n, prm, dx);	/* caller knows how to calc the current negative gradient, - df(x)/dxi  */
      esl_vec_DScale(dx, n, -1.0);
    } 
  else numeric_derivative(cfg, x, n, func, prm, dx, dat); /* else resort to brute force */
  esl_vec_DCopy(dx, n, cg);	/* and make that the first conjugate direction, cg  */


  /* (failsafe) convergence test: a completely zero direction can happen, 
   * and it either means we're stuck or we're finished (most likely stuck)
   */
  for (i1 = 0; i1 < n; i1++) 
    if (cg[i1] != 0.) break;
  if  (i1 == n) {
    if (opt_fx) *opt_fx = oldfx;
    free(wrk);
    return eslOK;
  }
  
  for (i = 1; i <= max_iterations; i++)
    {
      if (dat) {
	dat->niter    = i;  // this is how bracket(), brent(), and numeric_derivative() know what CG iteration they're on
	dat->nfunc[i] = 0;
      }
      
#if (eslDEBUGLEVEL >= 2)   // When debugging, it's useful to compare caller's deriv to numeric_deriv
      int j;
      printf("\nCG iteration %d\n", i+1);
      printf(" current point:       ");
      for (j = 0; j < n; j++) printf("%10.4g ", x[j]);
      printf("\n gradient:            ");
      for (j = 0; j < n; j++) printf("%10.4g ", dx[j]);
      numeric_derivative(cfg, x, n, func, prm, w1, dat);
      printf("\n numeric gradient:    ");
      for (j = 0; j < n; j++) printf("%10.4g ", w1[j]);
      printf("\n conjugate direction: ");
      for (j = 0; j < n; j++) printf("%10.4g ", cg[j]);
      printf("\n");
#endif

      /* Figure out the initial step size, <bx>, passed to bracketer. */
      if (u) bx = fabs(u[0] / cg[0]);
      else   bx = fabs(1.   / cg[0]);
      for (i1 = 1; i1 < n; i1++)
	{
	  if (u) cx = fabs(u[i1] / cg[i1]);
	  else   cx = fabs(1.    / cg[i1]);
	  if (cx < bx) bx = cx;
	}
 
      /* Bracket the minimum.*/
      bracket(cfg, x, cg, n, bx, func, prm, w1, &ax, &bx, &cx, &fa, &fb, &fc, dat);
       
      /* Minimize along the line given by the conjugate gradient <cg> */
      brent(cfg, x, cg, n, func, prm, ax, cx, w2, NULL, &fx, dat);
      esl_vec_DCopy(w2, n, x);

      /* Bail out if the function is now +/-inf: this can happen if the caller
       * has screwed something up.
       */
      if (! isfinite(fx)) ESL_XEXCEPTION(eslERANGE, "minimum not finite");

      /* Find the negative gradient at that point (temporarily in w1) */
      if (dfunc != NULL) 
	{
	  (*dfunc)(x, n, prm, w1);
	  esl_vec_DScale(w1, n, -1.0);
	}
      else numeric_derivative(cfg, x, n, func, prm, w1, dat); /* resort to brute force */

      /* Calculate the Polak-Ribiere coefficient */
      for (coeff = 0., i1 = 0; i1 < n; i1++)
	coeff += (w1[i1] - dx[i1]) * w1[i1];
      coeff /= esl_vec_DDot(dx, dx, n);
      
      /* Calculate the next conjugate gradient direction in w2 */
      esl_vec_DCopy(w1, n, w2);
      esl_vec_DAddScaled(w2, cg, coeff, n);

      /* Finishing set up for next iteration: */
      esl_vec_DCopy(w1, n, dx);
      esl_vec_DCopy(w2, n, cg);

      /* Now: x is the current point; 
       *      fx is the function value at that point;
       *      dx is the current gradient at x;
       *      cg is the current conjugate gradient direction. 
       */
      if (dat)
	dat->fx[i] = fx;

      /* Main convergence test. */
      if (esl_DCompare(fx, oldfx, cg_rtol, cg_atol) == eslOK) break;

      /* Second (failsafe) convergence test: a zero direction can happen, 
       * and it either means we're stuck or we're finished (most likely stuck)
       */
      for (i1 = 0; i1 < n; i1++) 
	if (cg[i1] != 0.) break;
      if  (i1 == n) break;

      oldfx = fx;
    }

  free(wrk);
  if (opt_fx) *opt_fx = fx;
  return (i > max_iterations ? eslENOHALT: eslOK);

 ERROR:
  free(wrk);
  if (opt_fx) *opt_fx = eslINFINITY;
  return status;
}



/*****************************************************************
 * 2. ESL_MIN_CFG: optional configuration/customization
 *****************************************************************/


/* Function:  esl_min_cfg_Create()
 * Synopsis:  Create customizable configuration for the minimizer
 * Incept:    SRE, Fri 20 Jul 2018 [Benasque]
 * 
 * Args:      n  - number of parameters that'll be optimized
 */
ESL_MIN_CFG *
esl_min_cfg_Create(int n)
{
  ESL_MIN_CFG *cfg = NULL;
  int          status;

  ESL_ALLOC(cfg, sizeof(ESL_MIN_CFG));
  cfg->u = NULL;
  ESL_ALLOC(cfg->u, sizeof(double) * n);
  
  cfg->max_iterations = eslMIN_MAXITER;
  cfg->cg_rtol        = eslMIN_CG_RTOL;
  cfg->cg_atol        = eslMIN_CG_ATOL;
  cfg->brent_rtol     = eslMIN_BRENT_RTOL;
  cfg->brent_atol     = eslMIN_BRENT_ATOL;
  cfg->brack_maxiter  = eslMIN_BRACK_MAXITER;
  cfg->deriv_step     = eslMIN_DERIV_STEP;
  cfg->n              = n;
  esl_vec_DSet(cfg->u, n, eslMIN_BRACK_STEP);
  return cfg;

 ERROR:
  esl_min_cfg_Destroy(cfg);
  return NULL;
}

/* Function:  esl_min_cfg_Destroy()
 * Incept:    SRE, Fri 20 Jul 2018 [Benasque]
 */
void
esl_min_cfg_Destroy(ESL_MIN_CFG *cfg)
{
  if (cfg)
    {
      free(cfg->u);
      free(cfg);
    }
}

/*****************************************************************
 * 3. ESL_MIN_DAT, for optional data collection
 *****************************************************************/

/* Function:  esl_min_dat_Create()
 * Synopsis:  Create an ESL_MIN_DAT for collecting stats from a run
 * Incept:    SRE, Fri 20 Jul 2018
 *
 * Purpose:   Collects statistics on a CG minimizer run for each
 *            iteration 1..niter, plus the initial position 0.
 *
 * Args:      cfg : optional custom config, or NULL.
 */
ESL_MIN_DAT *
esl_min_dat_Create(ESL_MIN_CFG *cfg)
{
  int max_iterations = cfg ? cfg->max_iterations : eslMIN_MAXITER;
  ESL_MIN_DAT *dat = NULL;
  int          status;

  ESL_ALLOC(dat, sizeof(ESL_MIN_DAT));
  dat->fx       = NULL;
  dat->brack_n  = NULL;
  dat->brack_ax = dat->brack_bx = dat->brack_cx = NULL;
  dat->brack_fa = dat->brack_fb = dat->brack_fc = NULL;
  dat->brent_n  = NULL;
  dat->brent_x  = NULL;
  dat->nfunc    = NULL;

  ESL_ALLOC(dat->fx,       sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_n,  sizeof(int)    * (max_iterations + 1));
  ESL_ALLOC(dat->brack_ax, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_bx, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_cx, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_fa, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_fb, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brack_fc, sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->brent_n,  sizeof(int)    * (max_iterations + 1));
  ESL_ALLOC(dat->brent_x,  sizeof(double) * (max_iterations + 1));
  ESL_ALLOC(dat->nfunc,    sizeof(int)    * (max_iterations + 1));

  /* initialize boundary conditions, unused vals at iteration 0 (initial point);
   * only dat->fx[0], ->nfunc[0] have data
   */
  dat->brack_n[0]  = 0;
  dat->brack_ax[0] = dat->brack_bx[0] = dat->brack_cx[0] = 0.;
  dat->brack_fa[0] = dat->brack_fb[0] = dat->brack_fc[0] = 0.;
  dat->brent_n[0]  = 0;
  dat->brent_x[0]  = 0.;

  return dat;

 ERROR:
  esl_min_dat_Destroy(dat);
  return NULL;
}


/* Function:  esl_min_dat_Destroy()
 * Synopsis:  Free an ESL_MIN_DAT
 * Incept:    SRE, Fri 20 Jul 2018 [Benasque]
 */
void
esl_min_dat_Destroy(ESL_MIN_DAT *dat)
{
  if (dat)
    {
      free(dat->fx);
      free(dat->brack_n);
      free(dat->brack_ax);
      free(dat->brack_bx);
      free(dat->brack_cx);
      free(dat->brack_fa);
      free(dat->brack_fb);
      free(dat->brack_fc);
      free(dat->brent_n);
      free(dat->brent_x);
      free(dat->nfunc);
      free(dat);
    }
}


/* Function:  esl_min_dat_Dump()
 * Synopsis:  Dump a table of stats from a minimizer run.
 * Incept:    SRE, Fri 20 Jul 2018 [Benasque]
 */
int
esl_min_dat_Dump(FILE *fp, ESL_MIN_DAT *dat)
{
  int iter;

  esl_dataheader(fp, 6, "iter", 16, "fx",       16, "diff",
		 7, "brack_n",  16, "brack_ax", 16, "brack_bx", 16, "brack_cx", 
		 16,"brack_fa", 16, "brack_fb", 16, "brack_fc",
		 7, "brent_n",  16, "brent_x",   5, "nfunc",
		 0);

  for (iter = 0; iter <= dat->niter; iter++)
    {
      fprintf(fp, "%6d %16g %16g %7d %16g %16g %16g %16g %16g %16g %7d %16g %5d\n",
	      iter,
	      dat->fx[iter],
	      iter > 0 ? dat->fx[iter-1] - dat->fx[iter] : 0.,
	      dat->brack_n[iter],
	      dat->brack_ax[iter], dat->brack_bx[iter], dat->brack_cx[iter],
	      dat->brack_fa[iter], dat->brack_fb[iter], dat->brack_fc[iter],
	      dat->brent_n[iter],  dat->brent_x[iter],  dat->nfunc[iter]);
    }
  return eslOK;
}


/*****************************************************************
 * 4. Internal functions: numeric deriv, bracketing, 1D line min
 *****************************************************************/


/* Return the negative gradient at a point, determined numerically.
 */
static void
numeric_derivative(ESL_MIN_CFG *cfg, double *x, int n, 
		   double (*func)(double *, int, void*),
		   void *prm, double *dx, ESL_MIN_DAT *dat)
{
  double  relstep = cfg ? cfg->deriv_step : eslMIN_DERIV_STEP;
  double *u       = cfg ? cfg->u          : NULL;
  int    i;
  double delta;
  double f1, f2;
  double tmp;

  for (i = 0; i < n; i++)
    {
      if (u) delta = fabs(u[i] * relstep);
      else   delta = fabs(relstep);

      tmp = x[i]; 
      x[i] = tmp + delta;
      f1  = (*func)(x, n, prm);
      x[i] = tmp - delta;
      f2  = (*func)(x, n, prm);
      x[i] = tmp;

      dx[i] = (-0.5 * (f1-f2)) / delta;

      if (dat) dat->nfunc[dat->niter] += 2;
      ESL_DASSERT1((! isnan(dx[i])));
    }
}

/* bracket():
 * SRE, Wed Jul 27 11:43:32 2005 [St. Louis]
 *
 * Purpose:   Bracket a minimum. 
 *
 *            The minimization is quasi-one-dimensional, 
 *            starting from an initial <n>-dimension vector <ori>
 *            in the <n>-dimensional direction <d>.
 *            
 *            Caller passes a ptr to the objective function <*func()>,
 *            and a void pointer to any necessary conditional 
 *            parameters <prm>. The objective function will
 *            be evaluated at a point <x> by calling
 *            <(*func)(x, n, prm)>. The caller's function
 *            is responsible to casting <prm> to whatever it's
 *            supposed to be, which might be a ptr to a structure,
 *            for example; typically, for a parameter optimization
 *            problem, this holds the observed data.
 *            
 *            The routine works in scalar multipliers relative
 *            to origin <ori> and direction <d>; that is, a new <n>-dimensional
 *            point <b> is defined as <ori> + <bx><d>, for a scalar <bx>.
 *            
 *            The routine identifies a triplet <ax>, <bx>, <cx> such
 *            that $a < b < c$ and such that a minimum is known to
 *            exist in the $(a,b)$ interval because $f(b) < f(a),
 *            f(c)$. Also, the <a..b> and <b...c> intervals are in
 *            a golden ratio; the <b..c> interval is 1.618 times larger
 *            than <a..b>.
 *
 *            Since <d> is usually in the direction of the gradient,
 *            the points <ax>,<bx>,<cx> might be expected to be $\geq 0$;
 *            however, when <ori> is already close to the minimum, 
 *            it is often faster to bracket the minimum using
 *            a negative <ax>. The caller might then try to be "clever"
 *            and assume that the minimum is in the <bx..cx> interval
 *            when <ax> is negative, rather than the full <ax..cx>
 *            interval. That cleverness can fail, though, if <ori>
 *            is already in fact the minimum, because the line minimizer
 *            in brent() assumes a non-inclusive interval. Use
 *            <ax..cx> as the bracket.
 *            
 * Args:      cfg       - OPTIONAL configuration of tunable params, or NULL
 *            ori       - n-dimensional starting vector
 *            d         - n-dimensional direction to minimize along
 *            n         - # of dimensions
 *            firststep - bx is initialized to this scalar multiplier
 *            *func()   - objective function to minimize
 *            prm       - void * to any constant data that *func() needs
 *            wrk       - workspace: 1 allocated n-dimensional vector
 *            ret_ax    - RETURN:  ax < bx < cx scalar bracketing triplet
 *            ret_bx    - RETURN:    ...ax may be negative
 *            ret_cx    - RETURN:    
 *            ret_fa    - RETURN:  function evaluated at a,b,c
 *            ret_fb    - RETURN:    ... f(b) < f(a),f(c)
 *            ret_fc    - RETURN:
 *            opt_dat   - optRETURN: stats collection, if caller provided for it
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENORESULT> if it fails to converge.
 *
 * Xref:      STL9/130.
 */
static int
bracket(ESL_MIN_CFG *cfg, double *ori, double *d, int n, double firststep,
	double (*func)(double *, int, void *), void *prm, 
	double *wrk, 
	double *ret_ax, double *ret_bx, double *ret_cx,
	double *ret_fa, double *ret_fb, double *ret_fc, ESL_MIN_DAT *dat)
{
  int    max_iterations = cfg ? cfg->brack_maxiter : eslMIN_BRACK_MAXITER;
  double ax,bx,cx;		/* scalar multipliers */
  double fa,fb,fc;		/* f() evaluations at those points */
  double swapper;
  int    niter;
  
  /* Set and evaluate our first two points f(a) and f(b), which
   * are initially at 0.0 and <firststep>.
   */
  ax = 0.;  /* always start w/ ax at the origin, ax=0 */
  fa = (*func)(ori, n, prm);

  bx = firststep;
  esl_vec_DCopy(ori, n, wrk);
  esl_vec_DAddScaled(wrk, d, bx, n);
  fb = (*func)(wrk, n, prm);

  /* In principle, we usually know that the minimum m lies to the
   * right of a, m>=a, because d is likely to be a gradient.  You
   * might think we want 0 = a < b < c.  In practice, there's problems
   * with that. It's far easier to identify bad points (f(x) > f(a))
   * than to identify good points (f(x) < f(a)), because letting f(x)
   * blow up to infinity is fine as far as bracketing is concerned.
   * It can be almost as hard to identify a point b that f(b) < f(a)
   * as it is to find the minimum in the first place!
   * Counterintuitively, in cases where f(b)>f(a), it's better
   * to just swap the a,b labels and look for c on the wrong side
   * of a! This often works immediately, if f(a) was reasonably
   * close to the minimum and f(b) and f(c) are both terrible.
   */
  if (fb > fa)
    {
      swapper = ax; ax = bx; bx = swapper;
      swapper = fa; fa = fb; fb = swapper;
    }

  /* Make our first guess at c.
   * Remember, we don't know that b>a any more, and c might go negative.
   * We'll either have:      a..b...c with a=0;
   *                or:  c...b..a     with b=0.
   * In many cases, we'll immediately be done.
   */
  cx = bx + (bx-ax)*1.618;
  esl_vec_DCopy(ori, n, wrk);
  esl_vec_DAddScaled(wrk, d, cx, n);
  fc = (*func)(wrk, n, prm);
  
  /* We're not satisfied until fb < fa, fc; 
   * throughout the routine, we guarantee that fb < fa;
   * so we just check fc.
   */
  niter = 0;
  while (fc <= fb)
    {
      /* Slide over, discarding the a point; choose 
       * new c point even further away.
       */
      ax = bx; bx = cx;
      fa = fb; fb = fc;
      cx = bx+(bx-ax)*1.618;
      esl_vec_DCopy(ori, n, wrk);
      esl_vec_DAddScaled(wrk, d, cx, n);
      fc = (*func)(wrk, n, prm);

      /* This is a rare instance. We've reach the minimum
       * by trying to bracket it. Also check that not all
       * three points are the same.
       */
      if (ax != bx && bx != cx && fa == fb && fb == fc) break;

      niter++;
      if (niter > max_iterations) ESL_EXCEPTION(eslENORESULT, "Failed to bracket a minimum.");
    }

  /* We're about to return. Assure the caller that the points
   * are in order a < b < c, not the other way.
   */
  if (ax > cx)
    {
      swapper = ax; ax = cx; cx = swapper;
      swapper = fa; fa = fc; fc = swapper;
    }

  if (dat)
    {
      dat->brack_n[dat->niter]  = niter;
      dat->nfunc[dat->niter]   += niter + 3;
      dat->brack_ax[dat->niter] = ax;
      dat->brack_bx[dat->niter] = bx;
      dat->brack_cx[dat->niter] = cx;
      dat->brack_fa[dat->niter] = fa;
      dat->brack_fb[dat->niter] = fb;
      dat->brack_fc[dat->niter] = fc;
    }

  /* Optional verbosity (depending on compile-time eslDEBUGLEVEL) */
  ESL_DPRINTF2(("bracket(), in %2d iterations, brackets a minimum: a,b,c = %10.4g %10.4g %10.4g\n", niter, ax, bx, cx));
  ESL_DPRINTF2(("bracket():                                     fa,fb,fc = %10.4g %10.4g %10.4g\n",        fa, fb, fc));

  /* Return. */
  *ret_ax = ax;  *ret_bx = bx;  *ret_cx = cx;
  *ret_fa = fa;  *ret_fb = fb;  *ret_fc = fc;
  return eslOK;
}

/* brent():
 * SRE, Sun Jul 10 19:07:05 2005 [St. Louis]
 *
 * Purpose:   Quasi-one-dimensional minimization of a function <*func()>
 *            in <n>-dimensions, along vector <dir> starting from a
 *            point <ori>. Identifies a scalar $x$ that approximates
 *            the position of the minimum along this direction, in a
 *            given bracketing interval (<a,b>).  The minimum must
 *            have been bracketed by the caller in the <(a,b)>
 *            interval.  <a> is often 0, because we often start at the
 *            <ori>.
 *
 *            A quasi-1D scalar coordinate $x$ (such as <a> or <b>) is
 *            transformed to a point $\mathbf{p}$ in n-space as:
 *            $\mathbf{p} = \mathbf{\mbox{ori}} + x
 *            \mathbf{\mbox{dir}}$.
 *
 *            Any extra (fixed) data needed to calculate <func> can be
 *            passed through the void <prm> pointer.
 *
 *            <dat->brent_rtol> and <dat->brent_atol> define the
 *            relative convergence tolerance, $\mbox{rtol} |x| +
 *            atol$. <rtol> should not be less than the square root of
 *            the machine precision.  The <DBL_EPSILON> is 2.2e-16 on
 *            many machines with 64-bit doubles, so <rtol> is on the
 *            order of 1e-8 or more. <atol> is a yet smaller number, used
 *            to avoid nonconvergence in the pathological case $x=0$.
 *
 *            Upon convergence (which is guaranteed), returns <xvec>,
 *            the n-dimensional minimum. Optionally, will also return
 *            <ret_x>, the scalar <x> that resulted in that
 *            n-dimensional minimum, and <ret_fx>, the objective
 *            function <*func(x)> at the minimum.
 *
 *            This is an implementation of the R.P. Brent (1973)
 *            algorithm for one-dimensional minimization without
 *            derivatives (modified from Brent's ALGOL60 code). Uses a
 *            combination of bisection search and parabolic
 *            interpolation; should exhibit superlinear convergence in
 *            most functions.
 *
 *
 * Args:      cfg     - OPTIONAL customization of tunable params, or NULL
 *            ori     - n-vector at origin
 *            dir     - direction vector (gradient) we're following from ori
 *            n       - dimensionality of ori, dir, and xvec
 *            (*func) - ptr to caller's objective function
 *            prm     - ptr to any additional data (*func)() needs
 *            a,b     - minimum is bracketed on interval [a,b]
 *            xvec    - RETURN: minimum, as an n-vector (caller allocated)
 *            opt_x   - optRETURN: scalar multiplier that gave xvec
 *            opt_fx  - optRETURN: f(x)
 *            dat     - optRETURN: statistics collection, or NULL
 *
 * Returns:   (void)
 *
 * Reference: See [Brent73], Chapter 5. My version is derived directly
 *            from Brent's description and his ALGOL60 code. I've
 *            preserved his variable names as much as possible, to
 *            make the routine follow his published description
 *            closely. The Brent algorithm is also discussed in
 *            Numerical Recipes [Press88].
 */
static void
brent(ESL_MIN_CFG *cfg, double *ori, double *dir, int n,
      double (*func)(double *, int, void *), void *prm, double a, double b, 
      double *xvec, double *opt_x, double *opt_fx, ESL_MIN_DAT *dat)
{
  double eps = cfg ? cfg->brent_rtol : eslMIN_BRENT_RTOL;
  double t   = cfg ? cfg->brent_atol : eslMIN_BRENT_ATOL;
  double tol;                   /* tolerance = eps|x| + t */
  double w,x,v,u;               /* with [a,b]: Brent's six points     */
  double m;                     /* midpoint of current [a,b] interval */
  double fu,fv,fw,fx;           /* function evaluations */
  double p,q;                   /* numerator, denominator of parabolic interpolation */
  double r;
  double d,e;                   /* last, next-to-last values of p/q  */
  double c = 1. - (1./eslCONST_GOLD); /* Brent's c; 0.381966; golden ratio */
  int    niter;			/* number of iterations */

  x=v=w= a + c*(b-a);           /* initial guess of x by golden section */
  esl_vec_DCopy(ori, n, xvec);  /* build xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  fx=fv=fw = (*func)(xvec, n, prm);   /* initial function evaluation */

  d = e = 0.;
  niter = 0;
  while (1) /* algorithm is guaranteed to converge. no maxiter needed. */
    {
      m   = 0.5 * (a+b);
      tol = eps*fabs(x) + t;
      if (fabs(x-m) <= 2*tol - 0.5*(b-a)) break; /* convergence test. */
      niter++;

      p = q = r = 0.;
      if (fabs(e) > tol)
        { /* Compute parabolic interpolation, u = x + p/q */
          r = (x-w)*(fx-fv);
          q = (x-v)*(fx-fw);
          p = (x-v)*q - (x-w)*r;
          q = 2*(q-r);
          if (q > 0) { p = -p; } else {q = -q;}
          r = e;
          e=d;                  /* e is now the next-to-last p/q  */
        }

      if (fabs(p) < fabs(0.5*q*r) || p < q*(a-x) || p < q*(b-x))
        { /* Seems well-behaved? Use parabolic interpolation to compute new point u */
          d = p/q;              /* d remembers last p/q */
          u = x+d;              /* trial point, for now... */

          if (2.0*(u-a) < tol || 2.0*(b-u) < tol) /* don't evaluate func too close to a,b */
            d = (x < m)? tol : -tol;
        }
      else /* Badly behaved? Use golden section search to compute u. */
        {
          e = (x<m)? b-x : a-x;  /* e = largest interval */
          d = c*e;
        }

      /* Evaluate f(), but not too close to x.  */
      if      (fabs(d) >= tol) u = x+d;
      else if (d > 0)          u = x+tol;
      else                     u = x-tol;
      esl_vec_DCopy(ori, n, xvec);  /* build xvec from ori, dir, u */
      esl_vec_DAddScaled(xvec, dir, u, n);
      fu = (*func)(xvec, n, prm);   /* f(u) */

      /* Bookkeeping.  */
     if (fu <= fx)
        {
          if (u < x) b = x; else a = x;
          v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        }
      else
        {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x)
            { v = w; fv = fw; w = u; fw = fu; }
          else if (fu <= fv || v==x || v ==w)
            { v = u; fv = fu; }
        }
    }

  if (dat)
    {
      dat->brent_n[dat->niter] = niter;
      dat->brent_x[dat->niter] = x;
      dat->nfunc [dat->niter] += niter + 1;
    }

  /* Optional verbosity, depending on compile-time eslDEBUGLEVEL */
  ESL_DPRINTF2(("  brent(), in %d iterations, finds 1D minimum at %10.4g, f() = %10.4g\n", niter, x, fx));

  /* Return */
  esl_vec_DCopy(ori, n, xvec);  /* build final xvec from ori, dir, x */
  esl_vec_DAddScaled(xvec, dir, x, n);
  if (opt_x)  *opt_x  = x;
  if (opt_fx) *opt_fx = fx;
}

/*****************************************************************
 * 5. Unit tests
 *****************************************************************/
#ifdef eslMINIMIZER_TESTDRIVE

/* f(x) = \sum_i a_i (x_i - b_i)^2 
 *    with an analytic minimum at x = b for a_i >= 0
 */
static double
test_func(double *x, int n, void *prm)
{
  double *a  =  (double *) prm; 
  double *b  = ((double *) prm) + n;
  double  fx = 0.;
  int     i;
  for (i = 0; i < n; i++) fx += a[i] * (x[i] - b[i]) * (x[i] - b[i]);
  return fx;
}
/* df(x)/dx_i = 2 a_i (x_i - b_i) */
static void
test_dfunc(double *x, int n, void *prm, double *dx)
{
  double *a  =  (double *) prm; 
  double *b  = ((double *) prm) + n; 
  int     i;
  for (i = 0; i < n; i++)
    dx[i] = 2. * a[i] * (x[i] - b[i]);
}

static void
utest_simplefunc(ESL_RANDOMNESS *rng)
{
  char    msg[] = "esl_minimizer simplefunc test failed";
  int     n     = 1 + esl_rnd_Roll(rng, 10);   // 1..10
  double *prm   = malloc(sizeof(double) * 2 * n);
  double *a     = prm;
  double *b     = prm+n;
  double *x     = malloc(sizeof(double) * n);
  double  fx;
  int     i;

  /* choose fixed params and an initial x[] randomly */
  for (i = 0; i < n; i++)
    {
      a[i] = esl_random(rng) * 10.;        // [0,10)
      b[i] = esl_random(rng) * 20. - 10.;  // [-10,10)
      x[i] = esl_random(rng) * 20. - 10.;  
    }

  esl_min_ConjugateGradientDescent(NULL, x, n, &test_func, &test_dfunc, (void *) prm, &fx, NULL);

  for (i = 0; i < n; i++)
    if ( esl_DCompare( b[i], x[i], 1e-5, 1e-10) != eslOK) esl_fatal(msg);
  if (esl_DCompare(0., fx, 0., 1e-5) != eslOK) esl_fatal(msg);  

  free(prm);
  free(x);
}

#endif /*eslMINIMIZER_TESTDRIVE*/


/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef eslMINIMIZER_TESTDRIVE
#include <esl_config.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_minimizer.h"
#include "esl_random.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                          docgroup*/
  { "-h",     eslARG_NONE,   NULL, NULL, NULL,  NULL,  NULL, NULL, "show brief help summary",             0 },
  { "-s",     eslARG_INT,    "42", NULL, NULL,  NULL,  NULL, NULL, "set random number generator seed",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = esl_getopts_CreateDefaultApp(options, 0, argc, argv,
						       "test driver for minimizer",
						       "[-options]");
  ESL_RANDOMNESS *rng = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));

  esl_fprintf(stderr, "## %s\n", argv[0]);
  esl_fprintf(stderr, "#  rng seed = %" PRIu32 "\n", esl_randomness_GetSeed(rng));

  utest_simplefunc(rng);  // test can stochastically fail - use fixed RNG seed in production code

  esl_fprintf(stderr, "#  status = ok\n");
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(rng);
  return 0;
}
#endif /*eslMINIMIZER_TESTDRIVE*/


/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef eslMINIMIZER_EXAMPLE
/*::cexcerpt::minimizer_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslMINIMIZER_EXAMPLE esl_minimizer.c esl_vectorops.c easel.c -lm
 * run:     ./example 
 */
#include <stdio.h>

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_minimizer.h"

/* a simple multidimensional quadratic w/ a minimum at 0:
 *    $f(x) = a_1 x_1^2 + ... a_n x_n^2$
 */ 
static double
example_func(double *x, int n, void *prm)
{
  double *a;
  double  fx;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (fx = 0., i = 0; i < n; i++)
    fx += a[i] * x[i] * x[i];
  return fx;
}
/* gradient of the f(x): d/dx_i = 2 a_i x_i
 */
static void
example_dfunc(double *x, int n, void *prm, double *dx)
{
  double *a;
  int     i;

  a = (double *) prm;	/* cast the data vector */
  for (i = 0; i < n; i++)
    dx[i] = 2.0 * a[i] * x[i];
}
int
main(int argc, char **argv)
{
  ESL_MIN_DAT *dat = esl_min_dat_Create(NULL);
  int    n = 6;
  double a[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double x[6] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
  double fx;
  int    i;

  esl_min_ConjugateGradientDescent(NULL, x, n, 
				   &example_func, &example_dfunc, (void *) a, 
				   &fx, dat);

  printf("At minimum: f(x) = %g\n", fx);
  printf("vector x = ");
  for (i = 0; i < 6; i++) printf("%g  ", x[i]);
  printf("\n");

  esl_min_dat_Dump(stdout, dat);

  esl_min_dat_Destroy(dat);
  return 0;
}
/*::cexcerpt::minimizer_example::end::*/
#endif /*eslMINIMIZER_EXAMPLE*/
