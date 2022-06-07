/*
 * Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>
 *
 * This file is licensed under the terms of the standard MIT License;
 * see the file LICENSE for details.
 */
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define Lx(i,j) ((((i)*((i)+1))>>1) + (j))

/* work: n values */
static bool LDLt_factor(
	double* A,
	const int n,
	double* work
	)
{
	bool result = true;
	double* v = work;
	for (int j=0; j<n; ++j)
	{
		for (int i = 0; i < j; ++i)
			v[i] = A[Lx(j,i)] * A[Lx(i,i)];
		v[j] = A[Lx(j,j)];
		for (int k = 0; k < j; ++k)
			v[j] -= A[Lx(j,k)] * v[k];
		if (v[j] <= 0.)
			result = false;
		A[Lx(j,j)] = v[j];
		for (int i = j + 1; i < n; ++i)
		{
			double s = A[Lx(i,j)];
			for (int k = 0; k < j; ++k)
				s -= A[Lx(i,k)] * v[k];
			A[Lx(i,j)] = s / v[j];
		}
	}
	return result;
}

/* work: n values */
static void LDLt_solve(
	const double* LD,
	const int n,
	double* b,
	double* work
	)
{
	double* zy = work;
	/* Forward; solve L.z = b */
	for (int i=0; i<n; ++i)
	{
		double s = b[i];
		for (int j=0; j<i; ++j)
			s -= LD[Lx(i,j)] * zy[j];
		zy[i] = s;
	}
	/* backward; solve D.y = z, L^t.x = y */
	for (int i=n-1; i>=0; --i)
	{
		double s = zy[i] / LD[Lx(i,i)];
		for (int k=i+1; k<n; ++k)
			s -= LD[Lx(k,i)] * b[k];
		b[i] = s;
	}
}

static inline double dot_v(
	const double* a,
	const double* b,
	const unsigned int n)
{
	double d = 0.;
	for (unsigned int k=0; k<n; ++k)
		d += a[k]*b[k];
	return d;
}

static inline double norm_v(
	const double* v,
	const unsigned int n)
{
	return sqrt(dot_v(v,v,n));
}

/* work: n values */
static inline double phi(
	double(*f)(const double*, void*),
	void* user_ptr,
	const double* xk,
	const double* pk,
	const double alpha,
	const unsigned int n,
	double* work
	)
{
	double* x = work;
	for (unsigned int i=0; i<n; ++i)
		x[i] = xk[i] + pk[i]*alpha;
	return f(x, user_ptr);
}

/* work: 2n values */
static inline double dphi(
	void(*Df)(const double*, double*, void*),
	void* user_ptr,
	const double* xk,
	const double* pk,
	double alpha,
	unsigned int n,
	double* work
	)
{
	double* x = work;
	double* Dfx = work + n;
	for (unsigned int i=0; i<n; ++i)
		x[i] = xk[i] + pk[i]*alpha;
	Df(x, Dfx, user_ptr);
	return dot_v(Dfx, pk, n);
}

/* work: 2n values */
static double zoom(
	double(*f)(const double*, void*),
	void(*Df)(const double*, double*, void*),
	void* user_ptr,
	const double* xk,
	const double* pk,
	const unsigned int n,
	const double c1,
	const double c2,
	double phi_0,
	double dphi_0,
	double alpha_lo,
	double alpha_hi,
	double* work
	)
{
	double alpha, phi_alpha;
	for (;;)
	{
		if (alpha_lo + 10*2.22e-16 >= alpha_hi)
			return alpha_lo;

		//printf("...zoom(%g,%g)\n", alpha_lo, alpha_hi);

		alpha = (alpha_lo + alpha_hi)*.5;
		phi_alpha = phi(f,user_ptr,xk,pk,alpha,n,work);
		if (phi_alpha > phi_0 + c1*alpha*dphi_0
			|| phi_alpha >= phi(f,user_ptr,xk,pk,alpha_lo,n,work))
		{
			alpha_hi = alpha;
		}
		else
		{
			double dphi_alpha = dphi(Df,user_ptr,xk,pk,alpha,n,work);
			if (fabs(dphi_alpha) <= -c2*dphi_0)
				return alpha;
			if (dphi_alpha*(alpha_hi - alpha_lo) >= 0.)
				alpha_hi = alpha_lo;
			alpha_lo = alpha;
		}
	}
}

/* work: 2n doubles */
static double linesearch(
	double(*f)(const double*, void*),
	void(*Df)(const double*, double*, void*),
	void* user_ptr,
	const double* xk,
	const double* pk,
	const unsigned int n,
	double* work
	)
{
	double alpha_prev = 0.;
	double alpha = 1.;

	double c1 = 1e-4;
	double c2 = 0.9;

	double phi_0 = phi(f,user_ptr,xk,pk,0.,n,work);
	double dphi_0 = dphi(Df,user_ptr,xk,pk,0.,n,work);
	double phi_alpha_prev;
	double alpha_star = 0.;
	double dphi_alpha;

	for (unsigned int iter=0; iter<1000; ++iter)
	{
		double phi_alpha = phi(f,user_ptr,xk,pk,alpha,n,work);
		if (phi_alpha > phi_0 + c1*alpha*dphi_0
			|| (iter>0 && phi_alpha >= phi_alpha_prev))
		{
			alpha_star = zoom(f,Df,user_ptr,xk,pk,n,c1,c2,
				phi_0,dphi_0,alpha_prev,alpha,work);
			break;
		}

		dphi_alpha = dphi(Df,user_ptr,xk,pk,alpha,n,work);
		if (fabs(dphi_alpha) >= 0.)
		{
			alpha_star = alpha;
			break;
		}

		if (dphi_alpha >= 0.)
		{
			alpha_star = zoom(f,Df,user_ptr,xk,pk,n,c1,c2,
				phi_0,dphi_0,alpha,alpha_prev,work);
			break;
		}

		alpha_prev = alpha;
		phi_alpha_prev = phi_alpha;
	}

	return alpha_star;
}

bool bfgs(
	double(*f)(const double*, void*),
	void(*Df)(const double*, double*, void*),
	void* user_ptr,
	double* x,
	const unsigned int n,
	const unsigned int max_iter
	)
{
	static const double eps = 2.2204e-16;
	const unsigned int ld = (n*(n+1))/2;
	double* B = (double*)malloc(ld*sizeof(double));
	double* LD = (double*)malloc(ld*sizeof(double));

	double* Dfx = (double*)malloc(n*sizeof(double));
	double* s = (double*)malloc(n*sizeof(double));
	double* y = (double*)malloc(n*sizeof(double));

	double* work = (double*)malloc(2*n*sizeof(double));
	double* Bs = work; // n values
	double* p = (double*)malloc(n*sizeof(double));
	double alpha;
	unsigned int iter;
	double step_length;

	for (iter=0; iter<max_iter; ++iter)
	{
		Df(x,Dfx,user_ptr);
		if (iter > 0)
		{
			for (unsigned int j=0; j<n; ++j)
			{
				s[j] += x[j];
				y[j] += Dfx[j];
			}

			double ys = dot_v(y,s,n);
			if (iter == 1)
			{
				double d = ys/dot_v(y,y,n);
				for (unsigned int i=0; i<n; ++i)
				{
					for (unsigned int j=0; j<i; ++j)
						B[Lx(i,j)] = 0.;
					B[Lx(i,i)] = d;
				}
			}

			/* B = B - ((B.s).(s'.B))/(s'.B.s) + (y.y')/(y'.s) */
			double sBs = 0.;
			for (unsigned int i=0; i<n; ++i)
			{
				Bs[i] = 0.;
				for (unsigned int j=0; j<=i; ++j)
					Bs[i] += B[Lx(i,j)]*s[j];
				for (unsigned int j=i+1; j<n; ++j)
					Bs[i] += B[Lx(j,i)]*s[j];

				sBs += s[i]*Bs[i];
			}
			double r_sBs = 1. / sBs;
			double r_ys = 1. / ys;
			for (unsigned int i=0; i<n; ++i)
				for (unsigned int j=0; j<=i; ++j)
					B[Lx(i,j)] -= Bs[i]*Bs[j]*r_sBs - y[i]*y[j]*r_ys;
		}

		if (norm_v(Dfx,n) < 100 * eps)
		{
			//printf("%3i: small gradient: |Dfx|=%g\n",
			//	iter, norm_v(Dfx,n));
			break;
		}

		for (unsigned int i=0; i<n; ++i)
			p[i] = Dfx[i];
		if (iter > 0)
		{
			memcpy(LD, B, ld*sizeof(double));
			if (LDLt_factor(LD, (int)n, work))
				LDLt_solve(LD, (int)n, p, work);
			/* else: singular Hessian approx.; trying gradient descent */
		}
		for (unsigned int i=0; i<n; ++i)
			p[i] = -p[i];

		alpha = linesearch(f,Df,user_ptr,x,p,n,work);
		step_length = alpha * norm_v(p,n);
		if (step_length < 10. * eps)
		{
			//printf("%3i: detected a small step (%g); aborting\n",
			//	iter, step_length);
			break;
		}
		//printf("%3i: obj. function: %g, step length: %g\n",
		//	iter, f(x,user_ptr), step_length);

		for (unsigned int i=0; i<n; ++i)
		{
			s[i] = -x[i];
			y[i] = -Dfx[i];
			x[i] += alpha * p[i];
		}
	}

	bool result = true;
	if (iter == max_iter)
	{
		/* max_iter exceeded */
		result = false;
	}
	bool pd = true;
	for (unsigned int i=0; i<n; ++i)
		if (LD[Lx(i,i)] <= 0.) { pd = false; break; }

	if (result && iter > 0 && !pd)
	{
		/* Hessian not positive definite; probably not a proper minimum */
		result = false;
	}

	free(p);
	free(work);
	free(y);
	free(s);
	free(Dfx);
	free(LD);
	free(B);

	return result;
}

// vim: fenc=utf-8 noet:
