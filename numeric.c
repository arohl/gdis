/*
Copyright (C) 2003 by Sean David Fleming
 
sean@ivec.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

The GNU GPL can also be found at http://www.gnu.org
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "gdis.h"
#include "numeric.h"

#define T_STEP 0.001
#define T_STEP_INV 1.0/T_STEP

#define ITMAX 100
#define EPS 3.0e-7

gdouble *s_data;
gint s_size=0;
gint sqrt_method=0;

/**********/
/* timing */
/**********/
gulong mytimer(void)
{
gulong usec;
struct timeval tv;

#ifndef __WIN32
gettimeofday(&tv, NULL);
#else
/* TODO */
#endif

usec = tv.tv_usec + 1000000*tv.tv_sec;

return(usec);
}

/************************************/
/* initialize the trig lookup table */
/************************************/
#define DEBUG_TRIG 0
void init_trig(void)
{
gint i;
gdouble a,s;

/* alloc */
s_size = PI/T_STEP + 1;
s_data = (gdouble *) g_malloc(s_size * sizeof(gdouble));

#if DEBUG_TRIG
printf("Initializing table: %d points\n", s_size);
#endif

i=0;
for (a=0.0 ; a<=PI ; a+= T_STEP)
  {
  s = sin(a);
  *(s_data+i) = s;
  i++;
  }
g_assert(i == s_size);
}

/********************/
/* SINE replacement */
/********************/
gdouble tbl_sin(gdouble angle)
{
gint quad, idx, idx2;
gdouble s,a,s1,s2,rem;

/* enforce range [0, 2PI] */
ca_rad(&angle);

/* determine quadrant */
a = angle;
quad = 0;
while (a > PI/2.0)
  {
  a -= PI/2.0;
  quad++;
  }
/* setup angle and sign accordingly */
s = 1.0;
switch (quad)
  {
  case 1:
    a += PI/2.0;
    break;
  case 3:
    a += PI/2.0;
  case 2:
    s *= -1.0;
    break;
  }
/* init retrieval indices */
idx2 = idx = (gint) (T_STEP_INV * a);
if (idx2)
  idx2--;
else
  idx2++;

g_assert(idx < s_size);
g_assert(idx2 >= 0);

/* get the two values */
s1 = *(s_data+idx) * s;
s2 = *(s_data+idx2) * s;
/* weighted average */
rem = T_STEP_INV * a - (gint) (T_STEP_INV * a);
s = (rem*s2 + (1.0-rem)*s1);

return(s);
}

/**********************/
/* COSINE replacement */
/**********************/
gdouble tbl_cos(gdouble angle)
{
return(tbl_sin(angle + PI/2.0));
}

/*******************************************/
/* get angle made with x axis of 2D vector */
/*******************************************/
gdouble angle_x_compute(gdouble x, gdouble y)
{
gdouble angle;

if (x == 0.0)
  x = G_MINDOUBLE;
if (fabs(y) <= 0.0001)
  y = G_MINDOUBLE;

angle = atan(y/x);

if (y >= 0.0 && x < 0.0)
  angle += G_PI;
else if (y < 0.0 && x < 0.0)
  angle += G_PI;
else if (y < 0.0 && x > 0.0)
  angle += 2.0*G_PI;

return(angle);
}

/***************************/
/* constrain angle [0,2pi] */
/***************************/
void ca_rad(gdouble *angle)
{
gint m;

m = *angle/(2.0*PI);

*angle -= (gdouble) m*2.0*PI;

if (*angle < 0.0)
  *angle += 2.0*PI;
}

/***************************/
/* constrain angle [0,360] */
/***************************/
void ca_deg(gdouble *angle)
{
gint m;

m = *angle/360.0;

*angle -= (gdouble) m*360.0;

if (*angle < 0.0)
  *angle += 360.0;
}

/********************/
/* SQRT replacement */
/********************/
gdouble fast_sqrt(gdouble s)
{
gdouble ds, r;

if (sqrt_method)
  {
/* NB: applicable range [0.1, 1.0] */
  r = 0.188030699 + 1.48359853*s - 1.0979059*s*s + 0.430357353*s*s*s;
  do
    {
    ds = s - r*r;
    r += 0.5*ds/r;
    }
  while (ds > FRACTION_TOLERANCE);
  }
else
  return(sqrt(s));

return(r);
}

/**************************************************/
/* test if approx sqrt is faster than system sqrt */
/**************************************************/
void init_sqrt(void)
{
/* TODO - test if it is faster */
sqrt_method = 1;
}

/*******************/
/* numerical setup */
/*******************/
void init_math(void)
{
init_trig();
init_sqrt();
}

/***************************************************/
/* convert a floating point to the nearest integer */
/***************************************************/
gint nearest_int(gdouble x)
{
gint i;

if (x > 0.0)
  i = (x+0.5);
else
  i = (x-0.5);
return(i);
}

/********************/
/* rounding routine */
/********************/
gdouble decimal_round(gdouble x, gint dp)
{
gint i;
gdouble y, f;

f = pow(10.0, dp);
y = x * f;
i = nearest_int(y);
y = i;
y /= f;
return(y);
}

/**********/
/* gammln */
/**********/
gdouble gammln(gdouble xx)
{
gint j;
gdouble x,tmp,ser;
static gdouble cof[6]={76.18009173,-86.50532033,24.01409822,
                       -1.231739516,0.120858003e-2,-0.536382e-5};
x = xx-1.0;
tmp = x+5.5;
tmp -= (x+0.5)*log(tmp);
ser = 1.0;
for (j=0 ; j<=5 ; j++)
  {
  x += 1.0;
  ser += cof[j]/x;
  }
return -tmp+log(2.50662827465*ser);
}

/*******/
/* gcd */
/*******/
gint gcd(gint p, gint q)
{
if (q == 0)
  return(abs(p));
else
  return(gcd(q, p%q));
}

/*******/
/* gcf */
/*******/
void gcf(gdouble *gammcf, gdouble a, gdouble x, gdouble *gln)
{
gint n;
gdouble gold=0.0,g,fac=1.0,b1=1.0;
gdouble b0=0.0,anf,ana,an,a1,a0=1.0;

*gln=gammln(a);
a1=x;
for (n=1 ; n<=ITMAX ; n++)
  {
  an = (gdouble) n;
  ana = an-a;
  a0 = (a1+a0*ana)*fac;
  b0 = (b1+b0*ana)*fac;
  anf = an*fac;
  a1 = x*a0+anf*a1;
  b1 = x*b0+anf*b1;
  if (a1) 
    {
    fac = 1.0/a1;
    g = b1*fac;
    if (fabs((g-gold)/g) < EPS)
      {
      *gammcf = exp(-x+a*log(x)-(*gln))*g;
      return;
      }
    gold = g;
    }
  }
g_assert_not_reached();
}

/********/
/* gser */
/********/
void gser(gdouble *gamser, gdouble a, gdouble x, gdouble *gln)
{
gint n;
gdouble sum,del,ap;

*gln = gammln(a);
if (x <= 0.0)
  {
  if (x < 0.0)
    g_assert_not_reached();
  *gamser=0.0;
  return;
  }
else
  {
  ap = a;
  del = sum=1.0/a;
  for (n=1 ; n<=ITMAX ; n++)
    {
    ap += 1.0;
    del *= x/ap;
    sum += del;
    if (fabs(del) < fabs(sum)*EPS)
      {
      *gamser = sum*exp(-x+a*log(x)-(*gln));
      return;
      }
    }
  g_assert_not_reached();
  return;
  }
}

/*********/
/* gammp */
/*********/
gdouble gammp(gdouble a, gdouble x)
{
gdouble gamser,gammcf,gln;

g_assert(x >= 0.0);
g_assert(a > 0.0);

if (x < (a+1.0))
  {
  gser(&gamser,a,x,&gln);
  return gamser;
  }
else
  {
  gcf(&gammcf,a,x,&gln);
  return 1.0-gammcf;
  }
}

/*********/
/* gammq */
/*********/
gdouble gammq(gdouble a, gdouble x)
{
gdouble gamser,gammcf,gln;

g_assert(x >= 0.0);
g_assert(a > 0.0);

if (x < (a+1.0))
  {
  gser(&gamser,a,x,&gln);
  return 1.0-gamser;
  }
else
  {
  gcf(&gammcf,a,x,&gln);
  return gammcf;
  }
}

/******************/
/* error function */
/******************/
gdouble erf(gdouble x)
{
gdouble val;

val = gammp(0.5, x*x);

if (x < 0.0)
  val *= -1.0;

return(val);
}

/********************************/
/* complementary error function */
/********************************/
gdouble erfc(gdouble x)
{
gdouble val;

if (x < 0.0)
  val = 1.0 + gammp(0.5, x*x);
else
  val = gammq(0.5, x*x);

return(val);
}

/****************/
/* general sort */
/****************/
void sort(gint size, gdouble *array)
{
gint i, swap;
gdouble tmp;

swap=1;
while (swap)
  {
  swap=0;
  for (i=1 ; i<size ; i++)
    {
/* TODO - direction flag for ascending or descending */
    if (array[i-1] > array[i])
      {
/* swap elements in array */
      tmp = array[i-1];
      array[i-1] = array[i];
      array[i] = tmp;
/* elements were swapped */
      swap++;
      }
    }
  }
}

/***************/
/* general min */
/***************/
gdouble min(gint size, gdouble *x)
{
gint i;
gdouble val;

g_assert(size > 0);

val = x[0];

for (i=1; i<size ; i++)
  if (x[i] < val)
    val = x[i];

return(val);
}

/***************/
/* general max */
/***************/
gdouble max(gint size, gdouble *x)
{
gint i;
gdouble val;

g_assert(size > 0);

val = x[0];

for (i=1; i<size ; i++)
  if (x[i] > val)
    val = x[i];

return(val);
}

/************************************/
/* numerical recipes - spline setup */
/************************************/
void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
int i, k;
double p, qn, sig, un, *u;

/* allocate 1 extra double as some people insist on starting at 1 */
u = g_malloc(n*sizeof(double));

if (yp1 > 0.99e30)
  y2[1] = u[1] = 0.0;
else
  {
  y2[1] = -0.5;
  u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }

for (i=2 ; i<=n-1 ; i++)
  {
  sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
  p = sig*y2[i-1]+2.0;
  y2[i] = (sig-1.0)/p;
  u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
  u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

if (ypn > 0.99e30)
  qn = un = 0.0;
else
  {
  qn = 0.5;
  un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }

y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0);

for (k=n-1 ; k>=1 ; k--)
  y2[k] = y2[k]*y2[k+1]+u[k];

g_free(u);
}

/********************************************/
/* numerical recipes - spline interpolation */
/********************************************/
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
int klo, khi, k;
double h, b, a;

klo = 1;
khi = n;

while (khi-klo > 1)
  {
  k = (khi+klo) >> 1;
  if (xa[k] > x)
    khi = k;
  else
    klo = k;
  }

h = xa[khi]-xa[klo];
if (h == 0.0)
  printf("splint(): bad xa input.\n");

a = (xa[khi]-x)/h;
b = (x-xa[klo])/h;
    
*y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h) / 6.0;
}

/****************************************/
/* numerical recipes - Cooley-Tukey FFT */
/****************************************/
#define SWAP(a, b) tempr=(a) ; (a) = (b) ; (b) = tempr;
void fft(gdouble *x, gint nn, gint isign)
{
gint n, mmax, m, j, istep, i;
gdouble wtemp, wr, wpr, wpi, wi, theta;
gdouble tempr, tempi;
gdouble *data = --x;    /* silly fortran programmers */

n = nn << 1;
j = 1;

for (i=1 ; i<n ; i+=2)
  {
  if (j > i)
    {
    SWAP(data[j], data[i]);
    SWAP(data[j+1], data[i+1]);
    }
  m = n >> 1;
  while (m >= 2 && j > m)
    {
    j -= m;
    m >>= 1;
    }
  j += m;
  }

mmax = 2;
while (n > mmax)
  {
  istep = 2 * mmax;
  theta = 2.0 * PI / (isign * mmax);
  wtemp = sin(0.5*theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0;
  wi = 0.0;
  for (m=1 ; m<mmax ; m+=2 )
    {
    for (i=m ; i<=n ; i+=istep)
      {
      j = i + mmax;
      tempr =  wr * data[j] - wi * data[j+1];
      tempi =  wr * data[j+1] - wi * data[j];
      data[j] = data[i] - tempr;
      data[j+1] = data[i+1] - tempi;
      data[i] += tempr;
      data[i+1] += tempi;
      }
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
    }
  mmax = istep;
  }
}

