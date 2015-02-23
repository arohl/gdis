
/* machine.h.  Generated automatically by configure.  */
/* Any machine specific stuff goes here */
/* Add details necessary for your own installation here! */

/* This is for use with "configure" -- if you are not using configure
	then use machine.van for the "vanilla" version of machine.h */

/* Note special macros: ANSI_C (ANSI C syntax)
			SEGMENTED (segmented memory machine e.g. MS-DOS)
			MALLOCDECL (declared if malloc() etc have
					been declared) */

/* #undef const */

/* #undef MALLOCDECL */
#define NOT_SEGMENTED 1
/* #undef HAVE_COMPLEX_H */

#define STDC_HEADERS 1
#define HAVE_BCOPY 1
#define HAVE_BZERO 1
/* #undef WORDS_BIGENDIAN */
#define U_INT_DEF 1
#define VARARGS 1


/* for basic or larger versions */
#define COMPLEX 1
#define SPARSE 1

/* for loop unrolling */
/* #undef VUNROLL */
/* #undef MUNROLL */

/* for segmented memory */
#ifndef NOT_SEGMENTED
#define	SEGMENTED
#endif

/* any compiler should have this header */
/* if not, change it */
#include        <stdio.h>

/* Check for ANSI C memmove and memset */
#ifdef STDC_HEADERS

/* standard copy & zero functions */
#define	MEM_COPY(from,to,size)	memmove((to),(from),(size))
#define	MEM_ZERO(where,size)	memset((where),'\0',(size))

#ifndef ANSI_C
#define ANSI_C 1
#endif

#endif

/* standard headers */
#include	<stdlib.h>
#include	<stddef.h>
#include	<string.h>
#include	<float.h>


/* if have bcopy & bzero and no alternatives yet known, use them */
#ifdef HAVE_BCOPY
#ifndef MEM_COPY
/* nonstandard copy function */
#define	MEM_COPY(from,to,size)	bcopy((char *)(from),(char *)(to),(int)(size))
#endif
#endif

#ifdef HAVE_BZERO
#ifndef MEM_ZERO
/* nonstandard zero function */
#define	MEM_ZERO(where,size)	bzero((char *)(where),(int)(size))
#endif
#endif

/* if the system has complex.h */
#ifdef HAVE_COMPLEX_H
#include	<complex.h>
#endif

/* If prototypes are available & ANSI_C not yet defined, then define it,
	but don't include any header files as the proper ANSI C headers
        aren't here */
#define HAVE_PROTOTYPES 1
#ifdef HAVE_PROTOTYPES
#ifndef ANSI_C
#define ANSI_C  1
#endif
#endif

/* floating point precision */

/* you can choose single, double or long double (if available) precision */

#define FLOAT 		1
#define DOUBLE 		2
#define LONG_DOUBLE 	3

/* #undef REAL_FLT */
/* #undef REAL_DBL */

/* if nothing is defined, choose double precision */
#ifndef REAL_DBL
#ifndef REAL_FLT
#define REAL_DBL 1
#endif
#endif

/* single precision */
#ifdef REAL_FLT
#define  Real float
#define  LongReal float
#define REAL FLOAT
#define LONGREAL FLOAT
#endif

/* double precision */
#ifdef REAL_DBL
#define Real double
#define LongReal double
#define REAL DOUBLE
#define LONGREAL DOUBLE
#endif


/* machine epsilon or unit roundoff error */
/* This is correct on most IEEE Real precision systems */
#ifdef DBL_EPSILON
#if REAL == DOUBLE
#define	MACHEPS	DBL_EPSILON
#elif REAL == FLOAT
#define	MACHEPS	FLT_EPSILON
#elif REAL == LONGDOUBLE
#define MACHEPS LDBL_EPSILON
#endif
#endif

#define F_MACHEPS 1.19209e-07
#define D_MACHEPS 2.22045e-16

#ifndef MACHEPS
#if REAL == DOUBLE
#define	MACHEPS	D_MACHEPS
#elif REAL == FLOAT  
#define MACHEPS F_MACHEPS
#elif REAL == LONGDOUBLE
#define MACHEPS D_MACHEPS
#endif
#endif

/* #undef M_MACHEPS */

/********************
#ifdef DBL_EPSILON
#define	MACHEPS	DBL_EPSILON
#endif
#ifdef M_MACHEPS
#ifndef MACHEPS
#define MACHEPS	M_MACHEPS
#endif
#endif
********************/

#define	M_MAX_INT 2147483647
#ifdef	M_MAX_INT
#ifndef MAX_RAND
#define	MAX_RAND ((double)(M_MAX_INT))
#endif
#endif

/* for non-ANSI systems */
#ifndef HUGE_VAL
#define HUGE_VAL HUGE
#endif


extern	int	isatty(int);



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


/* err.h  28/09/1993 */

/*  RCS id: $Id: mesch_pak.h,v 1.6 2006/09/25 05:11:51 seanfleming Exp $  */


#ifndef ERRHEADER
#define ERRHEADER


#include        <setjmp.h>

/* Error recovery */

extern	jmp_buf	restart;


/* max. # of error lists */
#define ERR_LIST_MAX_LEN   10

/* main error functions */

extern	int ev_err(const char *,int,int,const char *,int);  /* main error handler */
extern	int set_err_flag(int flag);         /* for different ways of handling
                                                errors, returns old value */
extern  int count_errs(int true_false);     /* to avoid "too many errors" */
extern  int err_list_attach(int list_num, int list_len,
	       char **err_ptr,int warn);  /* for attaching a list of errors */
extern  int err_is_list_attached(int list_num);  /* checking if a list 
						    is attached */
extern  int err_list_free(int list_num);   /* freeing a list of errors */


/* error(E_TYPE,"myfunc") raises error type E_TYPE for function my_func() */
#define	error(err_num,fn_name)	ev_err(__FILE__,err_num,__LINE__,fn_name,0)

/* warning(WARN_TYPE,"myfunc") raises warning type WARN_TYPE for 
   function my_func() */
#define warning(err_num,fn_name) ev_err(__FILE__,err_num,__LINE__,fn_name,1) 


/* error flags */
#define	EF_EXIT		0	/* exit on error */
#define	EF_ABORT	1	/* abort (dump core) on error */
#define	EF_JUMP		2	/* jump on error */
#define	EF_SILENT	3	/* jump, but don't print message */
#define	ERREXIT()	set_err_flag(EF_EXIT)
#define	ERRABORT()	set_err_flag(EF_ABORT)
/* don't print message */
#define	SILENTERR()	if ( ! setjmp(restart) ) set_err_flag(EF_SILENT)
/* return here on error */
#define	ON_ERROR()	if ( ! setjmp(restart) ) set_err_flag(EF_JUMP)


/* error types */
#define	E_UNKNOWN	0
#define	E_SIZES		1
#define	E_BOUNDS	2
#define	E_MEM		3
#define	E_SING		4
#define	E_POSDEF	5
#define	E_FORMAT	6
#define	E_INPUT		7
#define	E_NULL		8
#define	E_SQUARE	9
#define	E_RANGE		10
#define	E_INSITU2	11
#define	E_INSITU	12
#define	E_ITER		13
#define	E_CONV		14
#define	E_START		15
#define	E_SIGNAL	16
#define	E_INTERN	17
#define	E_EOF		18
#define E_SHARED_VECS   19
#define E_NEG           20
#define E_OVERWRITE     21
#define E_BREAKDOWN     22

/* warning types */
#define WARN_UNKNOWN     	0
#define WARN_WRONG_TYPE 	1
#define WARN_NO_MARK		2
#define WARN_RES_LESS_0         3
#define WARN_SHARED_VEC		4


/* error catching macros */

/* execute err_part if error errnum is raised while executing ok_part */
#define	catch(errnum,ok_part,err_part)	\
	{	jmp_buf _save;	int _err_num, _old_flag; \
		_old_flag = set_err_flag(EF_SILENT); \
		MEM_COPY(restart,_save,sizeof(jmp_buf)); \
		if ( (_err_num=setjmp(restart)) == 0 ) \
		{	ok_part; \
			set_err_flag(_old_flag); \
			MEM_COPY(_save,restart,sizeof(jmp_buf));	} \
		else if ( _err_num == errnum ) \
		{	set_err_flag(_old_flag);  \
			MEM_COPY(_save,restart,sizeof(jmp_buf)); \
			err_part;	} \
		else {	set_err_flag(_old_flag); \
			MEM_COPY(_save,restart,sizeof(jmp_buf)); \
			error(_err_num,"catch"); \
		} \
	}


/* execute err_part if any error raised while executing ok_part */
#define	catchall(ok_part,err_part) \
	{	jmp_buf _save;	int _err_num, _old_flag; \
		_old_flag = set_err_flag(EF_SILENT); \
		MEM_COPY(restart,_save,sizeof(jmp_buf)); \
		if ( (_err_num=setjmp(restart)) == 0 ) \
		{	ok_part; \
			set_err_flag(_old_flag); \
			MEM_COPY(_save,restart,sizeof(jmp_buf));	} \
		else \
		{	set_err_flag(_old_flag);  \
			MEM_COPY(_save,restart,sizeof(jmp_buf)); \
			err_part;	} \
	}


/* print message if error raised while executing ok_part,
                then re-raise error to trace calls */
#define	tracecatch(ok_part,function) \
	{	jmp_buf _save;	int _err_num, _old_flag; \
		_old_flag = set_err_flag(EF_JUMP); \
		MEM_COPY(restart,_save,sizeof(jmp_buf)); \
		if ( (_err_num=setjmp(restart)) == 0 ) \
		{	ok_part; \
			set_err_flag(_old_flag); \
			MEM_COPY(_save,restart,sizeof(jmp_buf));	} \
		else \
		{	set_err_flag(_old_flag);  \
			MEM_COPY(_save,restart,sizeof(jmp_buf)); \
			error(_err_num,function);	} \
	}



#endif   /* ERRHEADER */


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
		Type definitions for general purpose maths package
*/

#ifndef	MATRIXH

/* RCS id: $Id: mesch_pak.h,v 1.6 2006/09/25 05:11:51 seanfleming Exp $ */

#define	MATRIXH	


/* unsigned integer type */
/************************************************************
#ifndef U_INT_DEF
typedef	unsigned int	u_int;
#define U_INT_DEF
#endif
************************************************************/

/* vector definition */
typedef	struct	{
		unsigned int	dim, max_dim;
		Real	*ve;
		} VEC;

/* matrix definition */
typedef	struct	{
		unsigned int	m, n;
		unsigned int	max_m, max_n, max_size;
		Real	**me,*base;	/* base is base of alloc'd mem */
		} MAT;

/* band matrix definition */
typedef struct {
               MAT   *mat;       /* matrix */
               int   lb,ub;    /* lower and upper bandwidth */
               } BAND;


/* permutation definition */
typedef	struct	{
		unsigned int	size, max_size, *pe;
		} PERM;

/* integer vector definition */
typedef struct	{
		unsigned int	dim, max_dim;
		int	*ive;
	        } IVEC;


#ifndef MALLOCDECL
extern	void	*malloc(size_t),
		*calloc(size_t,size_t),
		*realloc(void *,size_t);
#endif /* MALLOCDECL */


#ifdef THREADSAFE
#define	STATIC
#else
#define	STATIC	static
#endif /* THREADSAFE */

void	m_version( void );

/* allocate one object of given type */
#define	NEW(type)	((type *)calloc((size_t)1,(size_t)sizeof(type)))

/* allocate num objects of given type */
#define	NEW_A(num,type)	((type *)calloc((size_t)(num),(size_t)sizeof(type)))

 /* re-allocate arry to have num objects of the given type */
#define	RENEW(var,num,type) \
    ((var)=(type *)((var) ? \
		    realloc((char *)(var),(size_t)((num)*sizeof(type))) : \
		    calloc((size_t)(num),(size_t)sizeof(type))))

#define	MEMCOPY(from,to,n_items,type) \
 MEM_COPY((char *)(from),(char *)(to),(unsigned)(n_items)*sizeof(type))


/* type independent min and max operations */
#ifndef max
#define	max(a,b)	((a) > (b) ? (a) : (b))
#endif /* max */
#ifndef min
#define	min(a,b)	((a) > (b) ? (b) : (a))
#endif /* min */


#undef TRUE
#define	TRUE	1
#undef FALSE
#define	FALSE	0


/* for input routines */
#define MAXLINE 81


/* Dynamic memory allocation */

/* Should use M_FREE/V_FREE/PX_FREE in programs instead of m/v/px_free()
   as this is considerably safer -- also provides a simple type check ! */



/* get/resize vector to given dimension */
extern	VEC *v_get(int), *v_resize(VEC *,int);
/* get/resize matrix to be m x n */
extern	MAT *m_get(int,int), *m_resize(MAT *,int,int);
/* get/resize permutation to have the given size */
/* get/resize an integer vector to given dimension */
/* get/resize a band matrix to given dimension */

/* free (de-allocate) (band) matrices, vectors, permutations and 
   integer vectors */
extern	int m_free(MAT *), v_free(VEC *), px_free(PERM *);



/* MACROS */

/* macros that also check types and sets pointers to NULL */
#define	M_FREE(mat)	( m_free(mat),	(mat)=(MAT *)NULL )
#define V_FREE(vec)	( v_free(vec),	(vec)=(VEC *)NULL )
#define	PX_FREE(px)	( px_free(px),	(px)=(PERM *)NULL )
#define	IV_FREE(iv)	( iv_free(iv),	(iv)=(IVEC *)NULL )

#define MAXDIM  	10000001


/* Entry level access to data structures */
/* routines to check indexes */
#define	m_chk_idx(A,i,j)	((i)>=0 && (i)<(A)->m && (j)>=0 && (j)<=(A)->n)
#define	v_chk_idx(x,i)		((i)>=0 && (i)<(x)->dim)
#define	bd_chk_idx(A,i,j)	((i)>=max(0,(j)-(A)->ub) && \
		(j)>=max(0,(i)-(A)->lb) && (i)<(A)->mat->n && (j)<(A)->mat->n)

#define	m_entry(A,i,j)		m_get_val(A,i,j)
#define	v_entry(x,i)		v_get_val(x,i)
#define	bd_entry(A,i,j)		bd_get_val(A,i,j)
#ifdef DEBUG
#define	m_set_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] = (val) : (error(E_BOUNDS,"m_set_val"), 0.0))
#define	m_add_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] += (val) : (error(E_BOUNDS,"m_add_val"), 0.0))
#define	m_sub_val(A,i,j,val)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] -= (val) : (error(E_BOUNDS,"m_sub_val"), 0.0))
#define	m_get_val(A,i,j)	( m_chk_idx(A,i,j) ? \
	(A)->me[(i)][(j)] : (error(E_BOUNDS,"m_get_val"), 0.0))
#define	v_set_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] = (val) : \
	(error(E_BOUNDS,"v_set_val"), 0.0))
#define	v_add_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] += (val) : \
	(error(E_BOUNDS,"v_set_val"), 0.0))
#define	v_sub_val(x,i,val)	( v_chk_idx(x,i) ? (x)->ve[(i)] -= (val) : \
	(error(E_BOUNDS,"v_set_val"), 0.0))
#define	v_get_val(x,i)	( v_chk_idx(x,i) ? (x)->ve[(i)] : \
	(error(E_BOUNDS,"v_get_val"), 0.0))
#define	bd_set_val(A,i,j,val)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] = (val) : \
	(error(E_BOUNDS,"bd_set_val"), 0.0))
#define	bd_add_val(A,i,j,val)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] += (val) : \
	(error(E_BOUNDS,"bd_set_val"), 0.0))
#define	bd_get_val(A,i,j)	( bd_chk_idx(A,i,j) ? \
	(A)->mat->me[(A)->lb+(j)-(i)][(j)] : \
	(error(E_BOUNDS,"bd_get_val"), 0.0))
#else /* no DEBUG */
#define	m_set_val(A,i,j,val)	((A)->me[(i)][(j)] = (val))
#define	m_add_val(A,i,j,val)	((A)->me[(i)][(j)] += (val))
#define	m_sub_val(A,i,j,val)	((A)->me[(i)][(j)] -= (val))
#define	m_get_val(A,i,j)	((A)->me[(i)][(j)])
#define	v_set_val(x,i,val)	((x)->ve[(i)] = (val))
#define	v_add_val(x,i,val)	((x)->ve[(i)] += (val))
#define	v_sub_val(x,i,val)	((x)->ve[(i)] -= (val))
#define	v_get_val(x,i)		((x)->ve[(i)])
#define	bd_set_val(A,i,j,val)	((A)->mat->me[(A)->lb+(j)-(i)][(j)] = (val))
#define	bd_add_val(A,i,j,val)	((A)->mat->me[(A)->lb+(j)-(i)][(j)] += (val))
#define	bd_get_val(A,i,j)	((A)->mat->me[(A)->lb+(j)-(i)][(j)])
#endif /* DEBUG */


/* I/O routines */

/* print x on file fp */
void v_foutput(FILE *fp,const VEC *x),
       /* print A on file fp */
	m_foutput(FILE *fp,const MAT *A),
       /* print px on file fp */
	px_foutput(FILE *fp,const PERM *px);
/* print ix on file fp */
void iv_foutput(FILE *fp,const IVEC *ix);

/* Note: if out is NULL, then returned object is newly allocated;
        Also: if out is not NULL, then that size is assumed */

/* read in vector from fp */
VEC *v_finput(FILE *fp,VEC *out);
/* read in matrix from fp */
MAT *m_finput(FILE *fp,MAT *out);
/* read in permutation from fp */
PERM *px_finput(FILE *fp,PERM *out);
/* read in int vector from fp */
IVEC *iv_finput(FILE *fp,IVEC *out);

/* fy_or_n -- yes-or-no to question in string s
        -- question written to stderr, input from fp 
        -- if fp is NOT a tty then return y_n_dflt */
int fy_or_n(FILE *fp, const char *s);

/* yn_dflt -- sets the value of y_n_dflt to val */
int yn_dflt(int val);

/* fin_int -- return integer read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, error exit otherwise
        -- ignore check if low > high           */
int fin_int(FILE *fp,const char *s,int low,int high);

/* fin_double -- return double read from file/stream fp
        -- prompt s on stderr if fp is a tty
        -- check that x lies between low and high: re-prompt if
                fp is a tty, error exit otherwise
        -- ignore check if low > high           */
double fin_double(FILE *fp,const char *s,double low,double high);

/* it skips white spaces and strings of the form #....\n
   Here .... is a comment string */
int skipjunk(FILE *fp);



/* MACROS */

/* macros to use stdout and stdin instead of explicit fp */
#define	v_output(vec)	v_foutput(stdout,vec)
#define	v_input(vec)	v_finput(stdin,vec)
#define	m_output(mat)	m_foutput(stdout,mat)
#define	m_input(mat)	m_finput(stdin,mat)
#define	px_output(px)	px_foutput(stdout,px)
#define	px_input(px)	px_finput(stdin,px)
#define	iv_output(iv)	iv_foutput(stdout,iv)
#define	iv_input(iv)	iv_finput(stdin,iv)

/* general purpose input routine; skips comments # ... \n */
#define	finput(fp,prompt,fmt,var) \
	( ( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) ), \
							fscanf(fp,fmt,var) )
#define	input(prompt,fmt,var)	finput(stdin,prompt,fmt,var)
#define	fprompter(fp,prompt) \
	( isatty(fileno(fp)) ? fprintf(stderr,prompt) : skipjunk(fp) )
#define	prompter(prompt)	fprompter(stdin,prompt)
#define	y_or_n(s)	fy_or_n(stdin,s)
#define	in_int(s,lo,hi)	fin_int(stdin,s,lo,hi)
#define	in_double(s,lo,hi)	fin_double(stdin,s,lo,hi)


/* special purpose access routines */

/* Copying routines */

/* copy in to out starting at out[i0][j0] */
extern	MAT	*_m_copy(const MAT *in,MAT *out,unsigned int i0,unsigned int j0),
		* m_move(const MAT *in, int, int, int, int, MAT *out, int, int),
		*vm_move(const VEC *in, int, MAT *out, int, int, int, int);
/* copy in to out starting at out[i0] */
extern	VEC	*_v_copy(const VEC *in,VEC *out,unsigned int i0),
		* v_move(const VEC *in, int, int, VEC *out, int),
		*mv_move(const MAT *in, int, int, int, int, VEC *out, int);
extern	PERM	*px_copy(const PERM *in,PERM *out);
extern	IVEC	*iv_copy(const IVEC *in,IVEC *out),
		*iv_move(const IVEC *in, int, int, IVEC *out, int);
extern  BAND    *bd_copy(const BAND *in,BAND *out);



/* MACROS */
#define	m_copy(in,out)	_m_copy(in,out,0,0)
#define	v_copy(in,out)	_v_copy(in,out,0)


/* Initialisation routines -- to be zero, ones, random or identity */
extern	VEC     *v_zero(VEC *), *v_rand(VEC *), *v_ones(VEC *);
extern	MAT     *m_zero(MAT *), *m_ident(MAT *), *m_rand(MAT *),
						*m_ones(MAT *);
extern	PERM    *px_ident(PERM *);
extern  IVEC    *iv_zero(IVEC *);

/* Basic vector operations */

extern	VEC	*sv_mlt(double s,const VEC *x,VEC *out),	/* out <- s.x */
		*mv_mlt(const MAT *A,const VEC *s,VEC *out),	/* out <- A.x */
		*vm_mlt(const MAT *A,const VEC *x,VEC *out),	/* out^T <- x^T.A */
		*v_add(const VEC *x,const VEC *y,VEC *out), 	/* out <- x + y */
                *v_sub(const VEC *x,const VEC *y,VEC *out),	/* out <- x - y */
		*px_vec(PERM *px,const VEC *x,VEC *out),	/* out <- P.x */
		*pxinv_vec(PERM *px,const VEC *x,VEC *out),	/* out <- P^{-1}.x */
		*v_mltadd(const VEC *x,const VEC *y,double s,VEC *out),   /* out <- x + s.y */
#ifdef PROTOTYPES_IN_STRUCT
		*v_map(double (*f)(double),const VEC *x,VEC *y),  
                                                 /* out[i] <- f(x[i]) */
		*_v_map(double (*f)(void *,double),void *p,const VEC *x,VEC *y),
#else
		*v_map(double (*f)(),const VEC *,VEC *), /* out[i] <- f(x[i]) */
		*_v_map(double (*f)(),void *,const VEC *,VEC *),
#endif /* PROTOTYPES_IN_STRUCT */
		*v_lincomb(int,const VEC **,const Real *,VEC *),   
                                                 /* out <- sum_i s[i].x[i] */
                *v_linlist(VEC *out,VEC *v1,double a1,...);
                                              /* out <- s1.x1 + s2.x2 + ... */

/* returns min_j x[j] (== x[i]) */
extern	double	v_min(const VEC *, int *), 
     /* returns max_j x[j] (== x[i]) */		
        v_max(const VEC *, int *), 
        /* returns sum_i x[i] */
        v_sum(const VEC *);

/* Hadamard product: out[i] <- x[i].y[i] */
extern	VEC	*v_star(const VEC *, const VEC *, VEC *),
                 /* out[i] <- x[i] / y[i] */
		*v_slash(const VEC *, const VEC *, VEC *),
               /* sorts x, and sets order so that sorted x[i] = x[order[i]] */ 
		*v_sort(VEC *, PERM *);

/* returns inner product starting at component i0 */
extern	double	_in_prod(const VEC *x, const VEC *y,unsigned int i0),
                /* returns sum_{i=0}^{len-1} x[i].y[i] */
                __ip__(const Real *,const Real *,int);

/* see v_mltadd(), v_add(), v_sub() and v_zero() */
extern	void	__mltadd__(Real *,const Real *,double,int),
		__add__(const Real *,const Real *,Real *,int),
		__sub__(const Real *,const Real *,Real *,int),
                __smlt__(const Real *,double,Real *,int),
		__zero__(Real *,int);


/* MACRO */
/* usual way of computing the inner product */
#define	in_prod(a,b)	_in_prod(a,b,0)

/* Norms */
/* scaled vector norms -- scale == NULL implies unscaled */
               /* returns sum_i |x[i]/scale[i]| */
extern	double	_v_norm1(const VEC *x,const VEC *scale),   
               /* returns (scaled) Euclidean norm */
                _v_norm2(const VEC *x,const VEC *scale),
               /* returns max_i |x[i]/scale[i]| */
		_v_norm_inf(const VEC *x,const VEC *scale);

/* unscaled matrix norms */
extern double m_norm1(const MAT *A), 
	m_norm_inf(const MAT *A), 
	m_norm_frob(const MAT *A);


/* MACROS */
/* unscaled vector norms */
#define	v_norm1(x)	_v_norm1(x,VNULL)
#define	v_norm2(x)	_v_norm2(x,VNULL)
#define	v_norm_inf(x)	_v_norm_inf(x,VNULL)

/* Basic matrix operations */

extern	MAT	*sm_mlt(double s, const MAT *A,MAT *out), 	/* out <- s.A */
		*m_mlt(const MAT *A,const MAT *B,MAT *out),	/* out <- A.B */
		*mmtr_mlt(const MAT *A,const MAT *B,MAT *out),	/* out <- A.B^T */
		*mtrm_mlt(const MAT *A,const MAT *B,MAT *out),	/* out <- A^T.B */
		*m_add(const MAT *A,const MAT *B,MAT *out),	/* out <- A + B */
		*m_sub(const MAT *A,const MAT *B,MAT *out),	/* out <- A - B */
		*sub_mat(const MAT *A,unsigned int,unsigned int,unsigned int,
			 unsigned int,MAT *out),
		*m_transp(const MAT *A,MAT *out),		/* out <- A^T */
                /* out <- A + s.B */ 
		*ms_mltadd(const MAT *A,const MAT *B,double s,MAT *out);   


extern  BAND    *bd_transp(const BAND *in, BAND *out),	/* out <- A^T */
  *sbd_mlt(Real s, const BAND *A, BAND *OUT),		/* OUT <- s.A */
  *bds_mltadd(const BAND *A, const BAND *B,double alpha, BAND *OUT),
  /* OUT <- A+alpha.B */
  *bd_zero(BAND *A);					/* A <- 0 */

extern	MAT	*px_rows(const PERM *px,const MAT *A,MAT *out),	/* out <- P.A */
		*px_cols(const PERM *px,const MAT *A,MAT *out),	/* out <- A.P^T */
		*swap_rows(MAT *,int,int,int,int),
		*swap_cols(MAT *,int,int,int,int),
                 /* A[i][j] <- out[j], j >= j0 */
		*_set_col(MAT *A,unsigned int i,const VEC *col,unsigned int j0),
                 /* A[i][j] <- out[i], i >= i0 */
		*_set_row(MAT *A,unsigned int j,const VEC *row,unsigned int i0);

extern	VEC	*get_row(const MAT *,unsigned int,VEC *),
		*get_col(const MAT *,unsigned int,VEC *),
		*sub_vec(const VEC *,int,int,VEC *),
                   /* mv_mltadd: out <- x + s.A.y */
		*mv_mltadd(const VEC *x,const VEC *y,const MAT *A,
			   double s,VEC *out),
                  /* vm_mltadd: out^T <- x^T + s.y^T.A */
		*vm_mltadd(const VEC *x,const VEC *y,const MAT *A,
			   double s,VEC *out),
                  /* bdv_mltadd: out <- x + s.A.y */
                *bdv_mltadd(const VEC *x,const VEC *y,const BAND *A,
			    double s,VEC *out);


/* MACROS */
/* row i of A <- vec */
#define	set_row(mat,row,vec)	_set_row(mat,row,vec,0) 
/* col j of A <- vec */
#define	set_col(mat,col,vec)	_set_col(mat,col,vec,0)


/* Basic permutation operations */

extern	PERM	*px_mlt(const PERM *px1,const PERM *px2,PERM *out),	/* out <- px1.px2 */
		*px_inv(const PERM *px,PERM *out),	/* out <- px^{-1} */
                 /* swap px[i] and px[j] */
		*px_transp(PERM *px,unsigned int i,unsigned int j);

     /* returns sign(px) = +1 if px product of even # transpositions
                           -1 if ps product of odd  # transpositions */
extern	int	px_sign(const PERM *);


/* Basic integer vector operations */

extern	IVEC	*iv_add(const IVEC *ix,const IVEC *iy,IVEC *out),  
  /* out <- ix + iy */
		*iv_sub(const IVEC *ix,const IVEC *iy,IVEC *out),  
  /* out <- ix - iy */
  /* sorts ix & sets order so that sorted ix[i] = old ix[order[i]] */
		*iv_sort(IVEC *ix, PERM *order);


/* miscellaneous functions */

double	square(double x), 	/* returns x^2 */
  cube(double x), 		/* returns x^3 */
  mrand(void);                  /* returns random # in [0,1) */

void	smrand(int seed),            /* seeds mrand() */
  mrandlist(Real *x, int len);       /* generates len random numbers */

void    m_dump(FILE *fp,const MAT *a), px_dump(FILE *fp, const PERM *px),
        v_dump(FILE *fp,const VEC *x), iv_dump(FILE *fp, const IVEC *ix);

MAT *band2mat(const BAND *bA, MAT *A);
BAND *mat2band(const MAT *A, int lb,int ub, BAND *bA);



/* miscellaneous constants */
#define	VNULL	((VEC *)NULL)
#define	MNULL	((MAT *)NULL)
#define	PNULL	((PERM *)NULL)
#define	IVNULL	((IVEC *)NULL)
#define BDNULL  ((BAND *)NULL)



/* varying number of arguments */

#include <stdarg.h>

/* prototypes */

int v_get_vars(int dim,...);
int iv_get_vars(int dim,...);
int m_get_vars(int m,int n,...);
int px_get_vars(int dim,...);

int v_resize_vars(int new_dim,...);
int iv_resize_vars(int new_dim,...);
int m_resize_vars(int m,int n,...);
int px_resize_vars(int new_dim,...);

int v_free_vars(VEC **,...);
int iv_free_vars(IVEC **,...);
int px_free_vars(PERM **,...);
int m_free_vars(MAT **,...);

#elif VARARGS
/* old varargs is used */

#include  <varargs.h>

/* prototypes */

int v_get_vars();
int iv_get_vars();
int m_get_vars();
int px_get_vars();

int v_resize_vars();
int iv_resize_vars();
int m_resize_vars();
int px_resize_vars();

int v_free_vars();
int iv_free_vars();
int px_free_vars();
int m_free_vars();


#endif /* MATRIXH */






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


/* meminfo.h  26/08/93 */
/* changed  11/12/93 */


#ifndef MEM_INFOH
#define MEM_INFOH



/* for hash table in mem_stat.c */
/* Note: the hash size should be a prime, or at very least odd */
#define MEM_HASHSIZE         509
#define MEM_HASHSIZE_FILE    "meminfo.h"


/* default: memory information is off */
/* set it to 1 if you want it all the time */
#define MEM_SWITCH_ON_DEF	0


/* available standard types */
#define TYPE_NULL              (-1)
#define TYPE_MAT    	        0
#define TYPE_BAND               1
#define TYPE_PERM		2
#define TYPE_VEC		3
#define TYPE_IVEC		4

#ifdef SPARSE
#define TYPE_ITER		5
#define TYPE_SPROW              6
#define TYPE_SPMAT		7
#endif

#ifdef COMPLEX
#ifdef SPARSE
#define TYPE_ZVEC		8
#define TYPE_ZMAT		9
#else
#define TYPE_ZVEC		5
#define TYPE_ZMAT		6
#endif
#endif

/* structure for memory information */
typedef struct {
   long bytes;       /* # of allocated bytes for each type (summary) */
   int  numvar;      /* # of allocated variables for each type */
} MEM_ARRAY;



#ifdef ANSI_C

int  mem_info_is_on(void);
int mem_info_on(int sw);

long mem_info_bytes(int type,int list);
int mem_info_numvar(int type,int list);
void mem_info_file(FILE * fp,int list);

void mem_bytes_list(int type,int old_size,int new_size,
		       int list);
void mem_numvar_list(int type, int num, int list);

#ifndef THREADSAFE
int mem_stat_reg_list(void **var,int type,int list,char *fname,int line);
int mem_stat_mark(int mark);
int mem_stat_free_list(int mark,int list);
int mem_stat_show_mark(void);
void mem_stat_dump(FILE *fp,int list);
int mem_attach_list(int list,int ntypes,char *type_names[],
	int (*free_funcs[])(), MEM_ARRAY info_sum[]);
int mem_free_vars(int list);
int mem_is_list_attached(int list);
void mem_dump_list(FILE *fp,int list);
int mem_stat_reg_vars(int list,int type,char *fname,int line,...);
#endif /* THREADSAFE */
#else
int mem_info_is_on();
int mem_info_on();

long mem_info_bytes();
int mem_info_numvar();
void mem_info_file();

void mem_bytes_list();
void mem_numvar_list();

#ifndef THREADSAFE
int mem_stat_reg_list();
int mem_stat_mark();
int mem_stat_free_list();
int mem_stat_show_mark();
void mem_stat_dump();
int mem_attach_list();
int mem_free_vars();
int mem_is_list_attached();
void mem_dump_list();
int mem_stat_reg_vars();
#endif /* THREADSAFE */

#endif 

/* macros */

#define mem_info()   mem_info_file(stdout,0)


#define mem_stat_reg(var,type)
#define MEM_STAT_REG(var,type)
#define mem_stat_free(mark)


/*
#ifndef THREADSAFE
#define mem_stat_reg(var,type)  mem_stat_reg_list((void **)var,type,0,__FILE__,__LINE__)
#define MEM_STAT_REG(var,type)  mem_stat_reg_list((void **)&(var),type,0,__FILE__,__LINE__)
#define mem_stat_free(mark)   mem_stat_free_list(mark,0)
#else
#define mem_stat_reg(var,type)
#define MEM_STAT_REG(var,type)
#define mem_stat_free(mark)
#endif
*/



#define mem_bytes(type,old_size,new_size)  \
  mem_bytes_list(type,old_size,new_size,0)

#define mem_numvar(type,num) mem_numvar_list(type,num,0)


/* internal type */

typedef struct {
   char **type_names;        /* array of names of types (strings) */
   int  (**free_funcs)();    /* array of functions for releasing types */
   unsigned ntypes;          /* max number of types */
   MEM_ARRAY *info_sum;      /* local array for keeping track of memory */
} MEM_CONNECT;

/* max number of lists of types */
#define MEM_CONNECT_MAX_LISTS    5


#endif

VEC	*symmeig(const MAT *, MAT *, VEC *);

