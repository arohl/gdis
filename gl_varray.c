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
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "numeric.h"
#include "opengl.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

/* externals */
extern struct sysenv_pak sysenv;


struct varray_pak
{
/* config */
gint gl_method;
gint periodic;
gint num_indices;
gint num_vertices;

/* data */
GLuint *indices;
/*
GLdouble *vertices;
GLdouble *new_vertices;
GLdouble *normals;
*/
GLfloat *interleaved;

GLdouble *colours;
};


/**********************/
/* creation primitive */
/**********************/
gpointer va_init(void)
{
struct varray_pak *va;

va = g_malloc(sizeof(struct varray_pak));

va->gl_method = 0;
va->periodic = FALSE;
va->num_indices = 0;
va->num_vertices = 0;

va->indices = NULL;
/*
va->vertices = NULL;
va->new_vertices = NULL;
va->normals = NULL;
*/
va->interleaved = NULL;

va->colours = NULL;

return(va);
}

/*************************/
/* destruction primitive */
/*************************/
void va_free(gpointer ptr_varray)
{
struct varray_pak *va = ptr_varray;

g_free(va->indices);
/*
g_free(va->vertices);
g_free(va->new_vertices);
g_free(va->normals);
*/

g_free(va->interleaved);

g_free(va->colours);
g_free(va);
}


/* FIXME - these should really be guints */
gint cyclic_clamp(gint n, gint start, gint stop)
{
gint range, offset;

range = stop - start + 1;
offset = n - start;
offset /= range;

return(n - offset*range);
}

gint cyclic_inc(gint n, gint start, gint stop)
{
if (n < stop)
  return(n+1);
else
  return(start);
}

gint cyclic_dec(gint n, gint start, gint stop)
{
if (n > start)
  return(n-1);
else
  return(stop);
}

/* CURRENT */
#define DEBUG_VA_GENERATE_INDICES 0
void va_generate_indices(gint depth, struct varray_pak *va)
{
gint d, i, j, m, n;
gint mdiv, div, mstart, nstart;

g_assert(va != NULL);

n = 18 * depth*depth;

#if DEBUG_VA_GENERATE_INDICES
printf("allocating for indices: %d\n", n);
#endif

va->num_indices = n;
va->indices = g_malloc(n * sizeof(GLuint));

i = 0;
n=1;
mstart = 0;
for (d=1 ; d<depth+1 ; d++)
  {
  nstart = n;
  div = d*6;

mdiv = 6*(d-1);

  m = mstart;

/* special case start */
  for (j=div ; j-- ; )
    {

*(va->indices+i++) = n;
#if DEBUG_VA_GENERATE_INDICES
printf(" %d", *(va->indices+i-1));
#endif

*(va->indices+i++) = cyclic_dec(n, nstart, nstart+div-1);
#if DEBUG_VA_GENERATE_INDICES
printf(" %d", *(va->indices+i-1));
#endif
n++;

*(va->indices+i++) = m++;
#if DEBUG_VA_GENERATE_INDICES
printf(" %d", *(va->indices+i-1));
#endif

    if ((j+1) % d == 0)
      m--;
    else
      {

/* missing triangle */
*(va->indices+i++) = cyclic_dec(m, mstart, nstart-1);
#if DEBUG_VA_GENERATE_INDICES
printf(" (%d", *(va->indices+i-1));
#endif
*(va->indices+i++) = cyclic_clamp(m, mstart, nstart-1);
#if DEBUG_VA_GENERATE_INDICES
printf(" %d", *(va->indices+i-1));
#endif
*(va->indices+i++) = cyclic_dec(n, nstart, nstart+div-1);
#if DEBUG_VA_GENERATE_INDICES
printf(" %d)", *(va->indices+i-1));
#endif
 
      }


    if (m >= nstart)
      {
      m -= nstart;
      m += mstart;
      }

    }

#if DEBUG_VA_GENERATE_INDICES
printf("\n");
#endif

  mstart = nstart;
  }

}

/*******************************/
/* sphere array initialization */
/*******************************/
#define DEBUG_MAKE_SPHERE 0
void va_make_sphere(gpointer ptr_varray)
{
gint i, j, k, n;
gdouble layers, divisions, alpha, beta, phi, theta, sp, cp, st, ct;
struct varray_pak *va = ptr_varray;

va->gl_method = GL_TRIANGLES;

/* number of layers - azimuthal sweep */
layers = sysenv.render.sphere_quality+1.0;
alpha = 0.5*PI/layers;

/* alloc */
n=0;
for (i=0 ; i<(gint) 1.0+layers ; i++)
  n += i;
n *= 6;
n++;

/* TODO - index & vertices depend on quality */
/*
sysenv.render.sphere_quality;
*/

/* setup indices */
va_generate_indices(layers, va);
/*
va->num_indices = 18;
va->indices = g_malloc(va->num_indices * sizeof(GLuint));
memcpy(va->indices, sphere_indices, va->num_indices*sizeof(GLuint));
*/

/* setup vertices and normals */
va->num_vertices = n;
/*
va->vertices = g_malloc(3 * n * sizeof(GLdouble));
va->new_vertices = g_malloc(3 * n * sizeof(GLdouble));
va->normals = g_malloc(3 * n * sizeof(GLdouble));
*/
va->interleaved = g_malloc(6 * n * sizeof(GLfloat));


/* initial vertex */
VEC3SET((va->interleaved), 0.0, 0.0, -1.0);
VEC3SET((va->interleaved+3), 0.0, 0.0, -1.0);

/* compute vertices */
k = 2;
for (i=1 ; i<(gint) 1.0+layers ; i++)
  {
#if DEBUG_MAKE_SPHERE
printf("[%d/%d]\n", i, (gint) layers);
#endif

  phi = i*alpha;
  sp = tbl_sin(phi);
  cp = tbl_cos(phi);

/* divisions at this level - hexagonal sweep */
  divisions = i*6.0;
  beta = 2.0*PI/divisions;
  for (j=0 ; j<(gint) divisions ; j++)
    {
    theta = j*beta;
    st = tbl_sin(theta);
    ct = tbl_cos(theta);

    VEC3SET((va->interleaved+3*k), ct*sp, st*sp, -cp);
    k++;
    VEC3SET((va->interleaved+3*k), ct*sp, st*sp, -cp);
    k++;

#if DEBUG_MAKE_SPHERE
printf("%7.4f %7.4f %7.4f\n", ct*sp, st*sp, -cp);
#endif
    }
  }

#if DEBUG_MAKE_SPHERE
printf("k = %d, n = %d\n", k, n);
#endif

/*
memcpy(va->normals, va->vertices, 3*n*sizeof(GLdouble));
*/
}

/*********************/
/* drawing primitive */
/*********************/
void va_draw_sphere(gpointer ptr_varray, gdouble *x, gdouble r)
{
gint i;
GLfloat *interleaved;
struct varray_pak *va = ptr_varray;

/* setup vertex/normal arrays */
interleaved = va->interleaved;
/*
memcpy(vertices, va->vertices, 3*va->num_vertices*sizeof(GLdouble));
*/

/* expand to required radius and move to desired location */
for (i=va->num_vertices ; i-- ; )
  {
  ARR3SET((interleaved+6*i), (va->interleaved+6*i+3));
  VEC3MUL((interleaved+6*i), r);
  ARR3ADD((interleaved+6*i), x);
  }

/* draw vertex array */
/*
glVertexPointer(3, GL_DOUBLE, 0, vertices);
if (va->normals)
  glNormalPointer(GL_DOUBLE, 0, va->normals);
*/

glVertexPointer(3, GL_FLOAT, 6*sizeof(GLfloat), &interleaved[0]);
glNormalPointer(GL_FLOAT, 6*sizeof(GLfloat), &interleaved[3]);


glDrawElements(va->gl_method, va->num_indices, GL_UNSIGNED_INT, va->indices);
}

