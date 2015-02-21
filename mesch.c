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

#include <stdio.h>
#include <glib.h>

#include "mesch_pak.h"

/*****************************/
/* matrix creation primitive */
/*****************************/
gpointer mesch_mat_new(gint i, gint j)
{
return(m_get(i,j));
}

/*****************************/
/* vector creation primitive */
/*****************************/
gpointer mesch_vec_new(gint i)
{
return(v_get(i));
}

/********************************/
/* matrix destruction primitive */
/********************************/
void mesch_m_free(gpointer data)
{
m_free(data);
}

/**********************/
/* element assignment */
/**********************/
void mesch_me_set(gpointer data, gint i, gint j, gdouble f)
{
MAT *mat = data;
mat->me[i][j] = f;
}

/********************/
/* element addition */
/********************/
void mesch_me_add(gpointer data, gint i, gint j, gdouble f)
{
MAT *mat = data;
mat->me[i][j] += f;
}

/**************************/
/* element multiplication */
/**************************/
void mesch_me_mul(gpointer data, gint i, gint j, gdouble f)
{
MAT *mat = data;
mat->me[i][j] *= f;
}

/****************************/
/* matrix element retrieval */
/****************************/
gdouble mesch_me_get(gpointer data, gint i, gint j)
{
MAT *mat = data;
return(mat->me[i][j]);
}

/****************************/
/* vector element retrieval */
/****************************/
gdouble mesch_ve_get(gpointer data, gint i)
{
VEC *vec = data;
return(vec->ve[i]);
}

/****************************/
/* vector element assignment */
/****************************/
void mesch_ve_set(gpointer data, gint i, gdouble value)
{
VEC *vec = data;
vec->ve[i] = value;
}

/***********************/
/* retrive matrix rows */
/***********************/
gint mesch_rows_get(gpointer data)
{
MAT *mat = data;
return(mat->m);
}

/***********************/
/* retrive matrix cols */
/***********************/
gint mesch_cols_get(gpointer data)
{
MAT *mat = data;
return(mat->n);
}

/***********************/
/* retrive vector dim  */
/***********************/
gint mesch_dim_get(gpointer data)
{
VEC *vec = data;
return(vec->dim);
}

/*******************************/
/* initialize a matrix to zero */
/*******************************/
void mesch_m_zero(gpointer data)
{
gint i, j;
MAT *mat = data;

for (i=mat->m ; i-- ; )
  for (j=mat->n ; j-- ; )
     mat->me[i][j] = 0.0;
}

/*******************************/
/* initialize a vector to zero */
/*******************************/
void mesch_v_zero(gpointer data)
{
gint i, n;

n = mesch_dim_get(data);
for (i=n ; i-- ; )
   mesch_ve_set(data, i, 0.0);
}

/************************************************************/
/* interface to eigenvalue calculation for symmetric matrix */
/************************************************************/
void mesch_sev_compute(gpointer m1, gpointer m2, gpointer v)
{
symmeig(m1, m2, v);
}
