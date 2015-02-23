
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

#include	<stdio.h>
#include	<math.h>
#include	<glib.h>

#include	"mesch_pak.h"


/* _m_copy -- copies matrix into new area
	-- out(i0:m,j0:n) <- in(i0:m,j0:n) */
MAT	*_m_copy(const MAT *in, MAT *out, unsigned int i0, unsigned int j0)
{
	unsigned int	i /* ,j */;

	if ( in==MNULL )
		g_assert_not_reached();
	if ( in==out )
		return (out);
	if ( out==MNULL || out->m < in->m || out->n < in->n )
		out = m_resize(out,in->m,in->n);

	for ( i=i0; i < in->m; i++ )
		MEM_COPY(&(in->me[i][j0]),&(out->me[i][j0]),
				(in->n - j0)*sizeof(Real));
		/* for ( j=j0; j < in->n; j++ )
			out->me[i][j] = in->me[i][j]; */

	return (out);
}

/* _v_copy -- copies vector into new area
	-- out(i0:dim) <- in(i0:dim) */
VEC	*_v_copy(const VEC *in, VEC *out, unsigned int i0)
{
	/* unsigned int	i,j; */

	if ( in==VNULL )
		g_assert_not_reached();
	if ( in==out )
		return (out);
	if ( out==VNULL || out->dim < in->dim )
		out = v_resize(out,in->dim);

	MEM_COPY(&(in->ve[i0]),&(out->ve[i0]),(in->dim - i0)*sizeof(Real));
	/* for ( i=i0; i < in->dim; i++ )
		out->ve[i] = in->ve[i]; */

	return (out);
}


/*
	The .._move() routines are for moving blocks of memory around
	within Meschach data structures and for re-arranging matrices,
	vectors etc.
*/

/* m_move -- copies selected pieces of a matrix
	-- moves the m0 x n0 submatrix with top-left cor-ordinates (i0,j0)
	   to the corresponding submatrix of out with top-left co-ordinates
	   (i1,j1)
	-- out is resized (& created) if necessary */
MAT	*m_move(const MAT *in, int i0,int j0, int m0,int n0,
		MAT *out, int i1, int j1)
{
    int		i;

    if ( ! in )
		g_assert_not_reached();
    if ( i0 < 0 || j0 < 0 || i1 < 0 || j1 < 0 || m0 < 0 || n0 < 0 ||
	 i0+m0 > in->m || j0+n0 > in->n )
		g_assert_not_reached();

    if ( ! out )
	out = m_resize(out,i1+m0,j1+n0);
    else if ( i1+m0 > out->m || j1+n0 > out->n )
	out = m_resize(out,max(out->m,i1+m0),max(out->n,j1+n0));

    for ( i = 0; i < m0; i++ )
	MEM_COPY(&(in->me[i0+i][j0]),&(out->me[i1+i][j1]),
		 n0*sizeof(Real));

    return out;
}

/* v_move -- copies selected pieces of a vector
	-- moves the length dim0 subvector with initial index i0
	   to the corresponding subvector of out with initial index i1
	-- out is resized if necessary */
VEC	*v_move(const VEC *in, int i0, int dim0,
		VEC *out, int i1)
{
    if ( ! in )
		g_assert_not_reached();
    if ( i0 < 0 || dim0 < 0 || i1 < 0 ||
	 i0+dim0 > in->dim )
		g_assert_not_reached();

    if ( (! out) || i1+dim0 > out->dim )
	out = v_resize(out,i1+dim0);

    MEM_COPY(&(in->ve[i0]),&(out->ve[i1]),dim0*sizeof(Real));

    return out;
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
		Files for matrix computations

	Givens operations file. Contains routines for calculating and
	applying givens rotations for/to vectors and also to matrices by
	row and by column.
*/

/* givens.c 1.2 11/25/87 */


/* givens -- returns c,s parameters for Givens rotation to
		eliminate y in the vector [ x y ]' */
void	givens(double x, double y, Real *c, Real *s)
{
	Real	norm;

	norm = sqrt(x*x+y*y);

	if ( norm == 0.0 )
	{	*c = 1.0;	*s = 0.0;	}	/* identity */
	else
	{	*c = x/norm;	*s = y/norm;	}
}

/* rot_vec -- apply Givens rotation to x's i & k components */
VEC	*rot_vec(const VEC *x,unsigned int i,unsigned int k, double c,double s,
		 VEC *out)
{
	Real	temp;

	if ( x==VNULL )
		g_assert_not_reached();
	if ( i >= x->dim || k >= x->dim )
		g_assert_not_reached();
	out = v_copy(x,out);

	/* temp = c*out->ve[i] + s*out->ve[k]; */
	temp = c*v_entry(out,i) + s*v_entry(out,k);
	/* out->ve[k] = -s*out->ve[i] + c*out->ve[k]; */
	v_set_val(out,k,-s*v_entry(out,i)+c*v_entry(out,k));
	/* out->ve[i] = temp; */
	v_set_val(out,i,temp);

	return (out);
}

/* rot_rows -- premultiply mat by givens rotation described by c,s */
MAT	*rot_rows(const MAT *mat, unsigned int i, unsigned int k,
		  double c, double s, MAT *out)
{
	unsigned int	j;
	Real	temp;

	if ( mat==(MAT *)NULL )
		g_assert_not_reached();
	if ( i >= mat->m || k >= mat->m )
		g_assert_not_reached();
	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->n; j++ )
	{
		/* temp = c*out->me[i][j] + s*out->me[k][j]; */
		temp = c*m_entry(out,i,j) + s*m_entry(out,k,j);
		/* out->me[k][j] = -s*out->me[i][j] + c*out->me[k][j]; */
		m_set_val(out,k,j, -s*m_entry(out,i,j) + c*m_entry(out,k,j));
		/* out->me[i][j] = temp; */
		m_set_val(out,i,j, temp);
	}

	return (out);
}

/* rot_cols -- postmultiply mat by givens rotation described by c,s */
MAT	*rot_cols(const MAT *mat,unsigned int i,unsigned int k,
		  gdouble c, gdouble s, MAT *out)
{
	unsigned int	j;
	Real	temp;

	if ( mat==(MAT *)NULL )
		g_assert_not_reached();
	if ( i >= mat->n || k >= mat->n )
		g_assert_not_reached();
	if ( mat != out )
		out = m_copy(mat,m_resize(out,mat->m,mat->n));

	for ( j=0; j<mat->m; j++ )
	{
/*
		 temp = c*out->me[j][i] + s*out->me[j][k]; 
		 out->me[j][k] = -s*out->me[j][i] + c*out->me[j][k]; 
		 out->me[j][i] = temp; 
*/
		temp = c*m_entry(out,j,i) + s*m_entry(out,j,k);
		m_set_val(out,j,k, -s*m_entry(out,j,i) + c*m_entry(out,j,k));
		m_set_val(out,j,i,temp);
	}

	return (out);
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/*
		Files for matrix computations

	Householder transformation file. Contains routines for calculating
	householder transformations, applying them to vectors and matrices
	by both row & column.
*/

/* hsehldr.c 1.3 10/8/87 */


/* hhvec -- calulates Householder vector to eliminate all entries after the
	i0 entry of the vector vec. It is returned as out. May be in-situ */
VEC	*hhvec(const VEC *vec, unsigned int i0, Real *beta,
	       VEC *out, Real *newval)
{
	Real	norm;

	out = _v_copy(vec,out,i0);
	norm = sqrt(_in_prod(out,out,i0));
	if ( norm <= 0.0 )
	{
		*beta = 0.0;
		return (out);
	}
	*beta = 1.0/(norm * (norm+fabs(out->ve[i0])));
	if ( out->ve[i0] > 0.0 )
		*newval = -norm;
	else
		*newval = norm;
	out->ve[i0] -= *newval;

	return (out);
}

/* hhtrvec -- apply Householder transformation to vector 
	-- that is, out <- (I-beta.hh(i0:n).hh(i0:n)^T).in
	-- may be in-situ */
VEC	*hhtrvec(const VEC *hh, double beta, unsigned int i0,
		 const VEC *in, VEC *out)
{
	Real	scale;
	/* unsigned int	i; */

	if ( hh==VNULL || in==VNULL )
		g_assert_not_reached();
	if ( in->dim != hh->dim )
		g_assert_not_reached();
	if ( i0 > in->dim )
		g_assert_not_reached();

	scale = beta*_in_prod(hh,in,i0);
	out = v_copy(in,out);
	__mltadd__(&(out->ve[i0]),&(hh->ve[i0]),-scale,(int)(in->dim-i0));
	/************************************************************
	for ( i=i0; i<in->dim; i++ )
		out->ve[i] = in->ve[i] - scale*hh->ve[i];
	************************************************************/

	return (out);
}

/* hhtrrows -- transform a matrix by a Householder vector by rows
	starting at row i0 from column j0 -- in-situ
	-- that is, M(i0:m,j0:n) <- M(i0:m,j0:n)(I-beta.hh(j0:n).hh(j0:n)^T) */
MAT	*hhtrrows(MAT *M, unsigned int i0, unsigned int j0,
		  const VEC *hh, double beta)
{
	Real	ip, scale;
	int	i /*, j */;

	if ( M==MNULL || hh==VNULL )
		g_assert_not_reached();
	if ( M->n != hh->dim )
		g_assert_not_reached();
	if ( i0 > M->m || j0 > M->n )
		g_assert_not_reached();

	if ( beta == 0.0 )	return (M);

	/* for each row ... */
	for ( i = i0; i < M->m; i++ )
	{	/* compute inner product */
		ip = __ip__(&(M->me[i][j0]),&(hh->ve[j0]),(int)(M->n-j0));
		/**************************************************
		ip = 0.0;
		for ( j = j0; j < M->n; j++ )
			ip += M->me[i][j]*hh->ve[j];
		**************************************************/
		scale = beta*ip;
		if ( scale == 0.0 )
		    continue;

		/* do operation */
		__mltadd__(&(M->me[i][j0]),&(hh->ve[j0]),-scale,
							(int)(M->n-j0));
		/**************************************************
		for ( j = j0; j < M->n; j++ )
			M->me[i][j] -= scale*hh->ve[j];
		**************************************************/
	}

	return (M);
}



/* _hhtrcols -- transform a matrix by a Householder vector by columns
	starting at row i0 from column j0 
	-- that is, M(i0:m,j0:n) <- (I-beta.hh(i0:m).hh(i0:m)^T)M(i0:m,j0:n)
	-- in-situ
	-- scratch vector w passed as argument
	-- raises error if w == NULL
*/
MAT	*_hhtrcols(MAT *M, unsigned int i0, unsigned int j0,
		   const VEC *hh, double beta, VEC *w)
{
	/* Real	ip, scale; */
	int	i /*, k */;
	/*  STATIC	VEC	*w = VNULL; */

	if ( M == MNULL || hh == VNULL || w == VNULL )
		g_assert_not_reached();
	if ( M->m != hh->dim )
		g_assert_not_reached();
	if ( i0 > M->m || j0 > M->n )
		g_assert_not_reached();

	if ( beta == 0.0 )	return (M);

	if ( w->dim < M->n )
	  w = v_resize(w,M->n);
	/*  MEM_STAT_REG(w,TYPE_VEC); */
	v_zero(w);

	for ( i = i0; i < M->m; i++ )
	    if ( hh->ve[i] != 0.0 )
		__mltadd__(&(w->ve[j0]),&(M->me[i][j0]),hh->ve[i],
							(int)(M->n-j0));
	for ( i = i0; i < M->m; i++ )
	    if ( hh->ve[i] != 0.0 )
		__mltadd__(&(M->me[i][j0]),&(w->ve[j0]),-beta*hh->ve[i],
							(int)(M->n-j0));
	return (M);
}



/* hhtrcols -- transform a matrix by a Householder vector by columns
	starting at row i0 from column j0 
	-- that is, M(i0:m,j0:n) <- (I-beta.hh(i0:m).hh(i0:m)^T)M(i0:m,j0:n)
	-- in-situ
	-- calls _hhtrcols() with the scratch vector w
	-- Meschach internal routines should call _hhtrcols() to
	avoid excessive memory allocation/de-allocation
*/
MAT	*hhtrcols(MAT *M, unsigned int i0, unsigned int j0,
		  const VEC *hh, double beta)
{
  STATIC VEC	*w = VNULL;

  if ( M == MNULL || hh == VNULL || w == VNULL )
		g_assert_not_reached();
  if ( M->m != hh->dim )
		g_assert_not_reached();
  if ( i0 > M->m || j0 > M->n )
		g_assert_not_reached();

  if ( ! w || w->dim < M->n )
    w = v_resize(w,M->n);
  MEM_STAT_REG(w,TYPE_VEC);

  M = _hhtrcols(M,i0,j0,hh,beta,w);

#ifdef THREADSAFE
  V_FREE(w);
#endif

  return M;
}


/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
		File containing routines for determining Hessenberg
	factorisations.
*/

/* Hfactor -- compute Hessenberg factorisation in compact form.
	-- factorisation performed in situ
	-- for details of the compact form see QRfactor.c and matrix2.doc */
MAT	*Hfactor(MAT *A, VEC *diag, VEC *beta)
{
	STATIC	VEC	*hh = VNULL, *w = VNULL;
	int	k, limit;

	if ( ! A || ! diag || ! beta )
		g_assert_not_reached();
	if ( diag->dim < A->m - 1 || beta->dim < A->m - 1 )
		g_assert_not_reached();
	if ( A->m != A->n )
		g_assert_not_reached();
	limit = A->m - 1;

	hh = v_resize(hh,A->m);
	w  = v_resize(w,A->n);
	MEM_STAT_REG(hh,TYPE_VEC);
	MEM_STAT_REG(w, TYPE_VEC);

	for ( k = 0; k < limit; k++ )
	  {
	    /* compute the Householder vector hh */
	    get_col(A,(unsigned int)k,hh);
	    /* printf("the %d'th column = ");	v_output(hh); */
	    hhvec(hh,k+1,&beta->ve[k],hh,&A->me[k+1][k]);
	    /* diag->ve[k] = hh->ve[k+1]; */
	    v_set_val(diag,k,v_entry(hh,k+1));
	    /* printf("H/h vector = ");	v_output(hh); */
	    /* printf("from the %d'th entry\n",k+1); */
	    /* printf("beta = %g\n",beta->ve[k]); */

	    /* apply Householder operation symmetrically to A */
	    _hhtrcols(A,k+1,k+1,hh,v_entry(beta,k),w);
	    hhtrrows(A,0  ,k+1,hh,v_entry(beta,k));
	    /* printf("A = ");		m_output(A); */
	  }

#ifdef THREADSAFE
	V_FREE(hh);	V_FREE(w);
#endif

	return (A);
}

/* makeHQ -- construct the Hessenberg orthogonalising matrix Q;
	-- i.e. Hess M = Q.M.Q'	*/
MAT	*makeHQ(MAT *H, VEC *diag, VEC *beta, MAT *Qout)
{
	int	i, j, limit;
	STATIC	VEC	*tmp1 = VNULL, *tmp2 = VNULL;

	if ( H==(MAT *)NULL || diag==(VEC *)NULL || beta==(VEC *)NULL )
		g_assert_not_reached();
	limit = H->m - 1;
	if ( diag->dim < limit || beta->dim < limit )
		g_assert_not_reached();
	if ( H->m != H->n )
		g_assert_not_reached();
	Qout = m_resize(Qout,H->m,H->m);

	tmp1 = v_resize(tmp1,H->m);
	tmp2 = v_resize(tmp2,H->m);
	MEM_STAT_REG(tmp1,TYPE_VEC);
	MEM_STAT_REG(tmp2,TYPE_VEC);

	for ( i = 0; i < H->m; i++ )
	{
		/* tmp1 = i'th basis vector */
		for ( j = 0; j < H->m; j++ )
			/* tmp1->ve[j] = 0.0; */
		    v_set_val(tmp1,j,0.0);
		/* tmp1->ve[i] = 1.0; */
		v_set_val(tmp1,i,1.0);

		/* apply H/h transforms in reverse order */
		for ( j = limit-1; j >= 0; j-- )
		{
			get_col(H,(unsigned int)j,tmp2);
			/* tmp2->ve[j+1] = diag->ve[j]; */
			v_set_val(tmp2,j+1,v_entry(diag,j));
			hhtrvec(tmp2,beta->ve[j],j+1,tmp1,tmp1);
		}

		/* insert into Qout */
		set_col(Qout,(unsigned int)i,tmp1);
	}

#ifdef THREADSAFE
	V_FREE(tmp1);	V_FREE(tmp2);
#endif

	return (Qout);
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
	This is a file of routines for zero-ing, and initialising
	vectors, matrices and permutations.
	This is to be included in the matrix.a library
*/


/* v_zero -- zero the vector x */
VEC	*v_zero(VEC *x)
{
	if ( x == VNULL )
		g_assert_not_reached();

	__zero__(x->ve,x->dim);
	/* for ( i = 0; i < x->dim; i++ )
		x->ve[i] = 0.0; */

	return x;
}

/* m_zero -- zero the matrix A */
MAT	*m_zero(MAT *A)
{
	int	i, A_m, A_n;
	Real	**A_me;

	if ( A == MNULL )
		g_assert_not_reached();

	A_m = A->m;	A_n = A->n;	A_me = A->me;
	for ( i = 0; i < A_m; i++ )
		__zero__(A_me[i],A_n);
		/* for ( j = 0; j < A_n; j++ )
			A_me[i][j] = 0.0; */

	return A;
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
  This file contains basic routines which are used by the functions
  in meschach.a etc.
  These are the routines that should be modified in order to take
  full advantage of specialised architectures (pipelining, vector
  processors etc).
  */


/* __ip__ -- inner product */
double	__ip__(const Real *dp1, const Real *dp2, int len)
{
#ifdef VUNROLL
    register int	len4;
    register Real	sum1, sum2, sum3;
#endif
    register int	i;
    register Real     sum;

    sum = 0.0;
#ifdef VUNROLL
    sum1 = sum2 = sum3 = 0.0;
    
    len4 = len / 4;
    len  = len % 4;
    
    for ( i = 0; i < len4; i++ )
    {
	sum  += dp1[4*i]*dp2[4*i];
	sum1 += dp1[4*i+1]*dp2[4*i+1];
	sum2 += dp1[4*i+2]*dp2[4*i+2];
	sum3 += dp1[4*i+3]*dp2[4*i+3];
    }
    sum  += sum1 + sum2 + sum3;
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	sum  += dp1[i]*dp2[i];
    
    return sum;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
void	__mltadd__(Real *dp1, const Real *dp2, double s, int len)
{
    register int	i;
#ifdef VUNROLL
    register int        len4;
    
    len4 = len / 4;
    len  = len % 4;
    for ( i = 0; i < len4; i++ )
    {
	dp1[4*i]   += s*dp2[4*i];
	dp1[4*i+1] += s*dp2[4*i+1];
	dp1[4*i+2] += s*dp2[4*i+2];
	dp1[4*i+3] += s*dp2[4*i+3];
    }
    dp1 += 4*len4;	dp2 += 4*len4;
#endif
    
    for ( i = 0; i < len; i++ )
	dp1[i] += s*dp2[i];
}

/* __smlt__ scalar multiply array c.f. sv_mlt() */
void	__smlt__(const Real *dp, double s, Real *out, int len)
{
    register int	i;
    for ( i = 0; i < len; i++ )
	out[i] = s*dp[i];
}

/* __add__ -- add arrays c.f. v_add() */
void	__add__(const Real *dp1, const Real *dp2, Real *out, int len)
{
    register int	i;
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] + dp2[i];
}

/* __sub__ -- subtract arrays c.f. v_sub() */
void	__sub__(const Real *dp1, const Real *dp2, Real *out, int len)
{
    register int	i;
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] - dp2[i];
}

/* __zero__ -- zeros an array of floating point numbers */
void	__zero__(Real *dp, int len)
{
gint i;

for ( i = len ; i-- ; )
  dp[i] = 0.0;
}


/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* memory.c 1.3 11/25/87 */


/* m_get -- gets an mxn matrix (in MAT form) by dynamic memory allocation
	-- normally ALL matrices should be obtained this way
	-- if either m or n is negative this will raise an error
	-- note that 0 x n and m x 0 matrices can be created */
MAT	*m_get(int m, int n)
{
   MAT	*matrix;
   int	i;
   
   if (m < 0 || n < 0)
		g_assert_not_reached();

   if ((matrix=NEW(MAT)) == (MAT *)NULL )
		g_assert_not_reached();
   
   matrix->m = m;		matrix->n = matrix->max_n = n;
   matrix->max_m = m;	matrix->max_size = m*n;

   if ((matrix->base = NEW_A(m*n,Real)) == (Real *)NULL )
   {
      free(matrix);
		g_assert_not_reached();
   }


   if ((matrix->me = (Real **)calloc(m,sizeof(Real *))) == 
       (Real **)NULL )
   {	free(matrix->base);	free(matrix);
		g_assert_not_reached();
     }
   
   /* set up pointers */
   for ( i=0; i<m; i++ )
     matrix->me[i] = &(matrix->base[i*n]);
   
   return (matrix);
}


/* v_get -- gets a VEC of dimension 'size'
   -- Note: initialized to zero */
VEC	*v_get(int size)
{
   VEC	*vector;
   
   if (size < 0)
		g_assert_not_reached();

   if ((vector=NEW(VEC)) == (VEC *)NULL )
		g_assert_not_reached();
   
   vector->dim = vector->max_dim = size;
   if ((vector->ve=NEW_A(size,Real)) == (Real *)NULL )
   {
      free(vector);
		g_assert_not_reached();
   }
   
   return (vector);
}

/* m_free -- returns MAT & asoociated memory back to memory heap */
int	m_free(MAT *mat)
{
 
   if ( mat==(MAT *)NULL || (int)(mat->m) < 0 ||
       (int)(mat->n) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( mat->base != (Real *)NULL ) {
      free((char *)(mat->base));
   }

   if ( mat->me != (Real **)NULL ) {
      free((char *)(mat->me));
   }
   
   free((char *)mat);
   
   return (0);
}



/* px_free -- returns PERM & asoociated memory back to memory heap */
int	px_free(PERM *px)
{
   if ( px==(PERM *)NULL || (int)(px->size) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( px->pe == (unsigned int *)NULL ) {
      free((char *)px);
   }
   else
   {
      free((char *)px->pe);
      free((char *)px);
   }
   
   return (0);
}



/* v_free -- returns VEC & asoociated memory back to memory heap */
int	v_free(VEC *vec)
{
   if ( vec==(VEC *)NULL || (int)(vec->dim) < 0 )
     /* don't trust it */
     return (-1);
   
   if ( vec->ve == (Real *)NULL ) {
      free((char *)vec);
   }
   else
   {
      free((char *)vec->ve);
      free((char *)vec);
   }
   
   return (0);
}



/* m_resize -- returns the matrix A of size new_m x new_n; A is zeroed
   -- if A == NULL on entry then the effect is equivalent to m_get() */
MAT	*m_resize(MAT *A,int new_m, int new_n)
{
   int	i;
   int	new_max_m, new_max_n, new_size, old_m, old_n;
   
   if (new_m < 0 || new_n < 0)
		g_assert_not_reached();

   if ( ! A )
     return m_get(new_m,new_n);

   /* nothing was changed */
   if (new_m == A->m && new_n == A->n)
     return A;

   old_m = A->m;	old_n = A->n;
   if ( new_m > A->max_m )
   {	/* re-allocate A->me */

      A->me = RENEW(A->me,new_m,Real *);
      if ( ! A->me )
		g_assert_not_reached();
   }
   new_max_m = max(new_m,A->max_m);
   new_max_n = max(new_n,A->max_n);
   
   new_size = new_max_m*new_max_n;
   if ( new_size > A->max_size )
   {	/* re-allocate A->base */

      A->base = RENEW(A->base,new_size,Real);
      if ( ! A->base )
		g_assert_not_reached();
      A->max_size = new_size;
   }
   
   /* now set up A->me[i] */
   for ( i = 0; i < new_m; i++ )
     A->me[i] = &(A->base[i*new_n]);
   
   /* now shift data in matrix */
   if ( old_n > new_n )
   {
      for ( i = 1; i < min(old_m,new_m); i++ )
	MEM_COPY((char *)&(A->base[i*old_n]),
		 (char *)&(A->base[i*new_n]),
		 sizeof(Real)*new_n);
   }
   else if ( old_n < new_n )
   {
      for ( i = (int)(min(old_m,new_m))-1; i > 0; i-- )
      {   /* copy & then zero extra space */
	 MEM_COPY((char *)&(A->base[i*old_n]),
		  (char *)&(A->base[i*new_n]),
		  sizeof(Real)*old_n);
	 __zero__(&(A->base[i*new_n+old_n]),(new_n-old_n));
      }
      __zero__(&(A->base[old_n]),(new_n-old_n));
      A->max_n = new_n;
   }
   /* zero out the new rows.. */
   for ( i = old_m; i < new_m; i++ )
     __zero__(&(A->base[i*new_n]),new_n);

   A->max_m = new_max_m;
   A->max_n = new_max_n;
   A->max_size = A->max_m*A->max_n;
   A->m = new_m;	A->n = new_n;
   
   return A;
}

/* v_resize -- returns the vector x with dim new_dim
   -- x is set to the zero vector */
VEC	*v_resize(VEC *x, int new_dim)
{
   
   if (new_dim < 0)
		g_assert_not_reached();

   if ( ! x )
     return v_get(new_dim);

   /* nothing is changed */
   if (new_dim == x->dim)
     return x;

   if ( x->max_dim == 0 )	/* assume that it's from sub_vec */
     return v_get(new_dim);
   
   if ( new_dim > x->max_dim )
   {

      x->ve = RENEW(x->ve,new_dim,Real);
      if ( ! x->ve )
		g_assert_not_reached();
      x->max_dim = new_dim;
   }
   
   if ( new_dim > x->dim )
     __zero__(&(x->ve[x->dim]),new_dim - x->dim);
   x->dim = new_dim;
   
   return x;
}


/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* 1.2 submat.c 11/25/87 */


/* get_col -- gets a specified column of a matrix and retruns it as a vector */
VEC	*get_col(const MAT *mat, unsigned int col, VEC *vec)
{
   unsigned int	i;
   
   if ( mat==(MAT *)NULL )
		g_assert_not_reached();
   if ( col >= mat->n )
		g_assert_not_reached();
   if ( vec==(VEC *)NULL || vec->dim<mat->m )
     vec = v_resize(vec,mat->m);
   
   for ( i=0; i<mat->m; i++ )
     vec->ve[i] = mat->me[i][col];
   
   return (vec);
}


/* _set_col -- sets column of matrix to values given in vec (in situ)
	-- that is, mat(i0:lim,col) <- vec(i0:lim) */
MAT	*_set_col(MAT *mat, unsigned int col, const VEC *vec, unsigned int i0)
{
   unsigned int	i,lim;
   
   if ( mat==(MAT *)NULL || vec==(VEC *)NULL )
		g_assert_not_reached();
   if ( col >= mat->n )
		g_assert_not_reached();
   lim = min(mat->m,vec->dim);
   for ( i=i0; i<lim; i++ )
     mat->me[i][col] = vec->ve[i];
   
   return (mat);
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
	File containing routines for symmetric eigenvalue problems
*/

#define	SQRT2	1.4142135623730949
#define	sgn(x)	( (x) >= 0 ? 1 : -1 )

/* trieig -- finds eigenvalues of symmetric tridiagonal matrices
	-- matrix represented by a pair of vectors a (diag entries)
		and b (sub- & super-diag entries)
	-- eigenvalues in a on return */
VEC	*trieig(VEC *a, VEC *b, MAT *Q)
{
	int	i, i_min, i_max, n, split;
	Real	*a_ve, *b_ve;
	Real	b_sqr, bk, ak1, bk1, ak2, bk2, z;
	Real	c, c2, cs, s, s2, d, mu;

	if ( ! a || ! b )
		g_assert_not_reached();
	if ( a->dim != b->dim + 1 || ( Q && Q->m != a->dim ) )
		g_assert_not_reached();
	if ( Q && Q->m != Q->n )
		g_assert_not_reached();

	n = a->dim;
	a_ve = a->ve;
	b_ve = b->ve;

	i_min = 0;
	while ( i_min < n )		/* outer while loop */
	{
		/* find i_max to suit;
			submatrix i_min..i_max should be irreducible */
		i_max = n-1;
		for ( i = i_min; i < n-1; i++ )
		    if ( b_ve[i] == 0.0 )
		    {	i_max = i;	break;	}
		if ( i_max <= i_min )
		{
		    /* printf("# i_min = %d, i_max = %d\n",i_min,i_max); */
		    i_min = i_max + 1;
		    continue;	/* outer while loop */
		}

		/* printf("# i_min = %d, i_max = %d\n",i_min,i_max); */

		/* repeatedly perform QR method until matrix splits */
		split = FALSE;
		while ( ! split )		/* inner while loop */
		{

		    /* find Wilkinson shift */
		    d = (a_ve[i_max-1] - a_ve[i_max])/2;
		    b_sqr = b_ve[i_max-1]*b_ve[i_max-1];
		    mu = a_ve[i_max] - b_sqr/(d + sgn(d)*sqrt(d*d+b_sqr));
		    /* printf("# Wilkinson shift = %g\n",mu); */

		    /* initial Givens' rotation */
		    givens(a_ve[i_min]-mu,b_ve[i_min],&c,&s);
		    s = -s;
		    /* printf("# c = %g, s = %g\n",c,s); */
		    if ( fabs(c) < SQRT2 )
		    {	c2 = c*c;	s2 = 1-c2;	}
		    else
		    {	s2 = s*s;	c2 = 1-s2;	}
		    cs = c*s;
		    ak1 = c2*a_ve[i_min]+s2*a_ve[i_min+1]-2*cs*b_ve[i_min];
		    bk1 = cs*(a_ve[i_min]-a_ve[i_min+1]) +
						(c2-s2)*b_ve[i_min];
		    ak2 = s2*a_ve[i_min]+c2*a_ve[i_min+1]+2*cs*b_ve[i_min];
		    bk2 = ( i_min < i_max-1 ) ? c*b_ve[i_min+1] : 0.0;
		    z  = ( i_min < i_max-1 ) ? -s*b_ve[i_min+1] : 0.0;
		    a_ve[i_min] = ak1;
		    a_ve[i_min+1] = ak2;
		    b_ve[i_min] = bk1;
		    if ( i_min < i_max-1 )
			b_ve[i_min+1] = bk2;
		    if ( Q )
			rot_cols(Q,i_min,i_min+1,c,-s,Q);
		    /* printf("# z = %g\n",z); */
		    /* printf("# a [temp1] =\n");	v_output(a); */
		    /* printf("# b [temp1] =\n");	v_output(b); */

		    for ( i = i_min+1; i < i_max; i++ )
		    {
			/* get Givens' rotation for sub-block -- k == i-1 */
			givens(b_ve[i-1],z,&c,&s);
			s = -s;

			/* perform Givens' rotation on sub-block */
		        if ( fabs(c) < SQRT2 )
		        {	c2 = c*c;	s2 = 1-c2;	}
		        else
		        {	s2 = s*s;	c2 = 1-s2;	}
		        cs = c*s;
			bk  = c*b_ve[i-1] - s*z;
			ak1 = c2*a_ve[i]+s2*a_ve[i+1]-2*cs*b_ve[i];
			bk1 = cs*(a_ve[i]-a_ve[i+1]) +
						(c2-s2)*b_ve[i];
			ak2 = s2*a_ve[i]+c2*a_ve[i+1]+2*cs*b_ve[i];
			bk2 = ( i+1 < i_max ) ? c*b_ve[i+1] : 0.0;

			z  = ( i+1 < i_max ) ? -s*b_ve[i+1] : 0.0;

			a_ve[i] = ak1;	a_ve[i+1] = ak2;
			b_ve[i] = bk1;
			if ( i < i_max-1 )
			    b_ve[i+1] = bk2;
			if ( i > i_min )
			    b_ve[i-1] = bk;
			if ( Q )
			    rot_cols(Q,i,i+1,c,-s,Q);
		        /* printf("# a [temp2] =\n");	v_output(a); */
		        /* printf("# b [temp2] =\n");	v_output(b); */
		    }

		    /* test to see if matrix should be split */
		    for ( i = i_min; i < i_max; i++ )
			if ( fabs(b_ve[i]) < MACHEPS*
					(fabs(a_ve[i])+fabs(a_ve[i+1])) )
			{   b_ve[i] = 0.0;	split = TRUE;	}

		    /* printf("# a =\n");	v_output(a); */
		    /* printf("# b =\n");	v_output(b); */
		}
	}

	return a;
}

/* symmeig -- computes eigenvalues of a dense symmetric matrix
	-- A **must** be symmetric on entry
	-- eigenvalues stored in out
	-- Q contains orthogonal matrix of eigenvectors
	-- returns vector of eigenvalues */
VEC	*symmeig(const MAT *A, MAT *Q, VEC *out)
{
	int	i;
	STATIC MAT	*tmp = MNULL;
	STATIC VEC	*b   = VNULL, *diag = VNULL, *beta = VNULL;

	if ( ! A )
		g_assert_not_reached();
	if ( A->m != A->n )
		g_assert_not_reached();
	if ( ! out || out->dim != A->m )
		out = v_resize(out,A->m);

	tmp  = m_resize(tmp,A->m,A->n);
	tmp  = m_copy(A,tmp);
	b    = v_resize(b,A->m - 1);
	diag = v_resize(diag,(unsigned int)A->m);
	beta = v_resize(beta,(unsigned int)A->m);
	MEM_STAT_REG(tmp,TYPE_MAT);
	MEM_STAT_REG(b,TYPE_VEC);
	MEM_STAT_REG(diag,TYPE_VEC);
	MEM_STAT_REG(beta,TYPE_VEC);

	Hfactor(tmp,diag,beta);
        if (Q)
          makeHQ(tmp,diag,beta,Q);

	for ( i = 0; i < A->m - 1; i++ )
	  {
		out->ve[i] = tmp->me[i][i];
		b->ve[i] = tmp->me[i][i+1];
	  }


	out->ve[i] = tmp->me[i][i];
	trieig(out,b,Q);

#ifdef	THREADSAFE
	M_FREE(tmp);	V_FREE(b);	V_FREE(diag);	V_FREE(beta);
#endif
	return out;
}

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/* vecop.c 1.3 8/18/87 */

/* _in_prod -- inner product of two vectors from i0 downwards
	-- that is, returns a(i0:dim)^T.b(i0:dim) */
double	_in_prod(const VEC *a, const VEC *b, unsigned int i0)
{
	unsigned int	limit;
	/* Real	*a_v, *b_v; */
	/* register Real	sum; */

	if ( a==(VEC *)NULL || b==(VEC *)NULL )
		g_assert_not_reached();
	limit = min(a->dim,b->dim);
	if ( i0 > limit )
		g_assert_not_reached();

	return __ip__(&(a->ve[i0]),&(b->ve[i0]),(int)(limit-i0));
	/*****************************************
	a_v = &(a->ve[i0]);		b_v = &(b->ve[i0]);
	for ( i=i0; i<limit; i++ )
		sum += a_v[i]*b_v[i];
		sum += (*a_v++)*(*b_v++);

	return (double)sum;
	******************************************/
}

