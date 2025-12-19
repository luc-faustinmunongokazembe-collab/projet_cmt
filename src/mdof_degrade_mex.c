/* mdof_degrade_mex.c  (corrected single-definition version)
 * Implements progressive-degradation MDOF solver:
 * - One-pass Newmark (beta=1/4, gamma=1/2 passed in)
 * - Rayleigh damping matrix C is constant and passed in
 * - Stiffness degrades based on instantaneous drift (damage permanent)
 * - Peak plastic drift tracked
 *
 * MATLAB signature:
 * [u,ud,udd,story_drifts,stiffness_history,deg_hist,plastic_hist,max_drift_ratio,yielded] = ...
 *   mdof_degrade_mex(M, C, k_elastic, story_heights, ag, dt, beta, gamma, ...
 *                    yield_drift, ultimate_drift, residual_strength, degradation_rate)
 */

#include "mex.h"
#include <math.h>
#include <string.h>

/* ---------------- helper: build shear building stiffness K from story k ---------------- */
static void build_shear_stiffness(const double *k_cur, int n, double *K)
{
    int i;
    /* zero */
    for (i = 0; i < n*n; ++i) K[i] = 0.0;
    /* assemble */
    for (i = 0; i < n; ++i) {
        K[i + i*n] += k_cur[i];
        if (i > 0) {
            K[i     + (i-1)*n] -= k_cur[i];
            K[(i-1) +  i    *n] -= k_cur[i];
            K[(i-1) + (i-1)*n] += k_cur[i];
        }
    }
}

/* ---------------- helper: small Ax=b solver with partial pivoting ---------------- */
static int solve_linear(double *A, double *b, int n)
{
    int i, j, k, piv;
    for (k = 0; k < n-1; ++k) {
        /* pivot on column k */
        piv = k;
        double maxabs = fabs(A[k + k*n]);
        for (i = k+1; i < n; ++i) {
            double v = fabs(A[i + k*n]);
            if (v > maxabs) { maxabs = v; piv = i; }
        }
        if (maxabs < 1e-20) return 1; /* singular */

        /* swap rows k <-> piv */
        if (piv != k) {
            for (j = k; j < n; ++j) {
                double tmp = A[k + j*n]; A[k + j*n] = A[piv + j*n]; A[piv + j*n] = tmp;
            }
            double tb = b[k]; b[k] = b[piv]; b[piv] = tb;
        }

        /* eliminate below */
        {
            double akk = A[k + k*n];
            for (i = k+1; i < n; ++i) {
                double factor = A[i + k*n] / akk;
                A[i + k*n] = 0.0;
                for (j = k+1; j < n; ++j)
                    A[i + j*n] -= factor * A[k + j*n];
                b[i] -= factor * b[k];
            }
        }
    }
    /* back substitution */
    for (i = n-1; i >= 0; --i) {
        double sum = b[i];
        for (j = i+1; j < n; ++j) sum -= A[i + j*n] * b[j];
        double aii = A[i + i*n];
        if (fabs(aii) < 1e-20) return 1;
        b[i] = sum / aii;
    }
    return 0;
}

/* ---------------- helper: compute story drifts from floor displacements ---------------- */
static void compute_drifts(const double *u, const double *h, int n, double *dr)
{
    int j;
    dr[0] = u[0] / h[0];
    for (j = 1; j < n; ++j)
        dr[j] = (u[j] - u[j-1]) / h[j];
}

/* --------------------- gateway --------------------- */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* ---- check I/O ---- */
    if (nrhs != 12) mexErrMsgTxt("Expected 12 inputs.");
    if (nlhs != 9)  mexErrMsgTxt("Expected 9 outputs.");

    /* inputs */
    const mxArray *M_in    = prhs[0];
    const mxArray *C_in    = prhs[1];
    const mxArray *kel_in  = prhs[2];
    const mxArray *h_in    = prhs[3];
    const mxArray *ag_in   = prhs[4];
    double dt              = mxGetScalar(prhs[5]);
    double beta            = mxGetScalar(prhs[6]);
    double gamma           = mxGetScalar(prhs[7]);
    double yield_drift     = mxGetScalar(prhs[8]);
    double ultimate_drift  = mxGetScalar(prhs[9]);
    double residual_strength = mxGetScalar(prhs[10]);
    double degradation_rate  = mxGetScalar(prhs[11]);

    int n  = (int)mxGetM(M_in);
    int nt = (int)mxGetNumberOfElements(ag_in);
    if ((int)mxGetN(M_in) != n || (int)mxGetM(C_in) != n || (int)mxGetN(C_in) != n)
        mexErrMsgTxt("M and C must be n-by-n.");
    if ((int)mxGetM(kel_in) != n || (int)mxGetN(kel_in) != 1)
        mexErrMsgTxt("k_elastic must be n-by-1.");
    if ((int)mxGetM(h_in) != n || (int)mxGetN(h_in) != 1)
        mexErrMsgTxt("story_heights must be n-by-1.");
    if (nt < 2) mexErrMsgTxt("ag must have at least 2 time samples.");

    const double *M   = mxGetPr(M_in);
    const double *C   = mxGetPr(C_in);
    const double *kel = mxGetPr(kel_in);
    const double *h   = mxGetPr(h_in);
    const double *agp = mxGetPr(ag_in);

    /* Newmark constants */
    double a0_nm = 1.0/(beta*dt*dt);
    double a1_nm = gamma/(beta*dt);
    double a2    = 1.0/(beta*dt);
    double a3    = 1.0/(2.0*beta) - 1.0;
    double a4    = gamma/beta - 1.0;
    double a5    = dt*(gamma/(2.0*beta) - 1.0);

    /* outputs */
    plhs[0] = mxCreateDoubleMatrix(n, nt, mxREAL); double *u   = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n, nt, mxREAL); double *ud  = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(n, nt, mxREAL); double *udd = mxGetPr(plhs[2]);

    plhs[3] = mxCreateDoubleMatrix(n, nt, mxREAL); double *story_drifts   = mxGetPr(plhs[3]);
    plhs[4] = mxCreateDoubleMatrix(n, nt, mxREAL); double *stiffness_hist = mxGetPr(plhs[4]);
    plhs[5] = mxCreateDoubleMatrix(n, nt, mxREAL); double *deg_hist       = mxGetPr(plhs[5]);
    plhs[6] = mxCreateDoubleMatrix(n, nt, mxREAL); double *plastic_hist   = mxGetPr(plhs[6]);
    plhs[7] = mxCreateDoubleMatrix(n, nt, mxREAL); double *maxdr_hist     = mxGetPr(plhs[7]);

    plhs[8] = mxCreateDoubleMatrix(n, 1,  mxREAL); double *yielded_out    = mxGetPr(plhs[8]);

    /* work arrays */
    double *k_cur = (double*)mxCalloc(n,   sizeof(double));
    double *deg_f = (double*)mxCalloc(n,   sizeof(double));
    double *plast = (double*)mxCalloc(n,   sizeof(double));
    double *dr    = (double*)mxCalloc(n,   sizeof(double));
    double *Ft    = (double*)mxCalloc(n,   sizeof(double));
    double *Peff  = (double*)mxCalloc(n,   sizeof(double));
    double *Kcur  = (double*)mxCalloc(n*n, sizeof(double));
    double *Keff  = (double*)mxCalloc(n*n, sizeof(double));
    double *rhs   = (double*)mxCalloc(n,   sizeof(double));

    /* initialize */
    int i, j;
    for (j = 0; j < n; ++j) {
        k_cur[j]       = kel[j];
        deg_f[j]       = 1.0;
        plast[j]       = 0.0;
        yielded_out[j] = 0.0;
    }
    /* histories at col 0; u/ud/udd already zero */
    for (j = 0; j < n; ++j) {
        stiffness_hist[j + 0*n] = kel[j];
        deg_hist[j + 0*n]       = 1.0;
        plastic_hist[j + 0*n]   = 0.0;
        maxdr_hist[j + 0*n]     = 0.0;
        story_drifts[j + 0*n]   = 0.0;
    }

    /* time stepping */
    for (i = 0; i < nt-1; ++i) {
        const double *u_i   = u  + i*n;
        const double *ud_i  = ud + i*n;
        const double *udd_i = udd + i*n;

        /* 1) update stiffness from instantaneous drifts of u(:,i) */
        compute_drifts(u_i, h, n, dr);
        for (j = 0; j < n; ++j) {
            double abs_d = dr[j] >= 0.0 ? dr[j] : -dr[j];
            if (abs_d > yield_drift) {
                yielded_out[j] = 1.0;

                /* --- track peak plastic drift --- */
                double incr = abs_d - yield_drift; if (incr < 0.0) incr = 0.0;
                if (incr > plast[j]) plast[j] = incr;

                /* --- exponential degradation --- */
                double ratio = abs_d / ultimate_drift; if (ratio > 1.0) ratio = 1.0;
                double g = exp(- degradation_rate * ratio * ratio);
                if (g < residual_strength) g = residual_strength;

                /* --- permanent damage --- */
                deg_f[j] = fmin(g, deg_f[j]);
                k_cur[j] = kel[j] * deg_f[j];
            } else {
                k_cur[j] = kel[j] * deg_f[j]; /* keep previous degraded stiffness */
            }
        }

        /* 2) assemble tangent and Keff */
        build_shear_stiffness(k_cur, n, Kcur);

        /* Keff = Kcur + a0*M + a1*C */
        for (j = 0; j < n*n; ++j) Keff[j] = Kcur[j];
        for (j = 0; j < n;   ++j) Keff[j + j*n] += a0_nm * M[j + j*n];
        for (j = 0; j < n*n; ++j) Keff[j] += a1_nm * C[j];

        /* 3) effective load */
        {
            double ag_next = agp[i+1];
            for (j = 0; j < n; ++j)
                Ft[j] = - M[j + j*n] * ag_next;
        }
        for (j = 0; j < n; ++j) {
            double Mu = M[j + j*n]*( a0_nm*u_i[j] + a2*ud_i[j] + a3*udd_i[j] );
            double Cu = 0.0;
            int k;
            for (k = 0; k < n; ++k)
                Cu += C[j + k*n]*( a1_nm*u_i[k] + a4*ud_i[k] + a5*udd_i[k] );
            Peff[j] = Ft[j] + Mu + Cu;
        }

        /* 4) solve Keff * u(:,i+1) = Peff */
        for (j = 0; j < n; ++j) rhs[j] = Peff[j];
        {
            double *A = (double*)mxCalloc(n*n, sizeof(double));
            memcpy(A, Keff, n*n*sizeof(double));
            if (solve_linear(A, rhs, n) != 0) {
                mxFree(A);
                mexErrMsgTxt("Linear solve failed (singular Keff).");
            }
            mxFree(A);
        }

        /* write u(:,i+1) */
        {
            double *u_ip1 = u + (i+1)*n;
            for (j = 0; j < n; ++j) u_ip1[j] = rhs[j];
        }

        /* 5) velocities & accelerations at i+1 (Newmark) */
        {
            const double *u_ip1 = u + (i+1)*n;
            double *ud_ip1  = ud  + (i+1)*n;
            double *udd_ip1 = udd + (i+1)*n;
            for (j = 0; j < n; ++j) {
                ud_ip1[j]  = a1_nm*(u_ip1[j] - u_i[j]) - a4*ud_i[j] - a5*udd_i[j];
                udd_ip1[j] = a0_nm*(u_ip1[j] - u_i[j]) - a2*ud_i[j] - a3*udd_i[j];
            }
        }

        /* 6) store histories at i+1 */
        for (j = 0; j < n; ++j) {
            story_drifts[j + (i+1)*n]   = dr[j];
            stiffness_hist[j + (i+1)*n] = k_cur[j];
            deg_hist[j + (i+1)*n]       = deg_f[j];
            plastic_hist[j + (i+1)*n]   = plast[j];

            /* max drift envelope */
            {
                double prev = maxdr_hist[j + i*n];
                double cur  = dr[j]; if (cur < 0.0) cur = -cur;
                maxdr_hist[j + (i+1)*n] = (prev > cur) ? prev : cur;
            }
        }
    } /* end time loop */

    /* free */
    mxFree(k_cur);
    mxFree(deg_f);
    mxFree(plast);
    mxFree(dr);
    mxFree(Ft);
    mxFree(Peff);
    mxFree(Kcur);
    mxFree(Keff);
    mxFree(rhs);
}