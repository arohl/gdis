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
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "morph.h"
#include "model.h"
#include "matrix.h"
#include "quaternion.h"
#include "render.h"
#include "space.h"
#include "surface.h"
#include "gui_shorts.h"
#include "interface.h"


#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/************************/
/* save morphology data */
/************************/
gint write_gmf(gchar *filename, struct model_pak *data)
{
GSList *plist, *slist;
struct plane_pak *pdata;
struct shift_pak *sdata;
FILE *fp;

g_return_val_if_fail(data != NULL, 1);
/* TODO - flag visible / option to save or not? */

/* open for saving */
fp = fopen(filename,"wt");
if (!fp)
  return(2);

/* header */
fprintf(fp, "\n title: ");
gdis_blurb(fp);
fprintf(fp, "  name: %s\n", data->basename);
fprintf(fp, " space: %s\n", data->sginfo.spacename);
fprintf(fp, "  cell: %f %f %f  %f %f %f\n",
             data->pbc[0], data->pbc[1], data->pbc[2],
             180.0*data->pbc[3]/PI, 180.0*data->pbc[4]/PI, 180.0*data->pbc[5]/PI);
/* default display type */
switch(data->morph_type)
  {
  case EQUIL_UN:
    fprintf(fp, " morph: unrelaxed equilibrium\n\n");
    break;
  case EQUIL_RE:
    fprintf(fp, " morph: relaxed equilibrium\n\n");
    break;
  case GROWTH_UN:
    fprintf(fp, " morph: unrelaxed growth\n\n");
    break;
  case GROWTH_RE:
    fprintf(fp, " morph: relaxed growth\n\n");
    break;
  default:
    fprintf(fp, " morph: Dhkl based\n\n");
    break;
  }

/* print out the planes */
for (plist=data->planes ; plist ; plist=g_slist_next(plist))
  {
  pdata = plist->data;
  if (!pdata->primary)
    continue;

  fprintf(fp,"miller:  %3d %3d %3d\n", 
              pdata->index[0], pdata->index[1], pdata->index[2]);

  for (slist=pdata->shifts ; slist ; slist=g_slist_next(slist)) 
    {
    sdata = (struct shift_pak *) slist->data;
  
    fprintf(fp, " %8.6f  %2d %2d  %10.4f %10.4f %10.4f %10.4f  %g\n",
            sdata->shift, sdata->region[0], sdata->region[1], 
            sdata->esurf[0], sdata->eatt[0],
            sdata->esurf[1], sdata->eatt[1], sdata->gnorm);
    }
  }

fclose(fp);
return(0);
}

/******************************************/
/* load miller & normal lengths from file */
/******************************************/
/* FIXME - when importing planes, this can mess up the periodic bond calc */
#define DEBUG_LOAD_PLANES 0
gint load_planes(gchar *filename, struct model_pak *data)
{
gint h, k, l, num_tokens;
gint cflag=FALSE, sflag=FALSE, gflag=FALSE;
gdouble m[3];
gchar **buff;
GSList *list, *new_planes;
struct plane_pak *plane=NULL;
struct shift_pak *shift=NULL;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* get next line */
new_planes = NULL;
for (;;)
  {
  buff = get_tokenized_line(fp, &num_tokens);
  if (!buff)
    break;

/* NB: only update space/cell etc. data if this call did */
/* not originate from an import planes call */
  if (data->id == MORPH)
    {
/* cell parameters */
    if (g_ascii_strncasecmp(*buff,"cell",4) == 0)
      {
      if (num_tokens >= 7)
        {
        cflag=TRUE;
        data->pbc[0] = str_to_float(*(buff+1));
        data->pbc[1] = str_to_float(*(buff+2));
        data->pbc[2] = str_to_float(*(buff+3));
        data->pbc[3] = PI*str_to_float(*(buff+4))/180.0;
        data->pbc[4] = PI*str_to_float(*(buff+5))/180.0;
        data->pbc[5] = PI*str_to_float(*(buff+6))/180.0;
/* compute direct & reciprocal lattices */
/* NB: enables fn_make_plane() to correctly compute Dhkl */
        matrix_lattice_init(data);
        }
      else
        printf("load_planes() error: bad cell line.\n");
      }

/* space group */
    if (g_ascii_strncasecmp(*buff,"space",5) == 0)
      {
      sflag=TRUE;
      data->sginfo.spacename = g_strjoinv(" ", buff+1);
      data->sginfo.spacenum = 0;
      }

/* default morphology type */
    if (g_ascii_strncasecmp(*buff, "morph", 5) == 0)
      {
      if (num_tokens >= 3)
        {
        if (g_ascii_strncasecmp(*(buff+1), "unrelaxed", 9) == 0)
          {
          if (g_ascii_strncasecmp(*(buff+2), "equil", 5) == 0)
            data->morph_type = EQUIL_UN;
          if (g_ascii_strncasecmp(*(buff+2), "growth", 6) == 0)
            data->morph_type = GROWTH_UN;
          }
  
        if (g_ascii_strncasecmp(*(buff+1), "relaxed", 7) == 0)
          {
          if (g_ascii_strncasecmp(*(buff+2), "equil", 5) == 0)
            data->morph_type = EQUIL_RE;
          if (g_ascii_strncasecmp(*(buff+2), "growth", 6) == 0)
            data->morph_type = GROWTH_RE;
          }
        }
      else
        printf("load_planes() error: bad type line.\n");
      }
    }

/* process miller line */
  if (g_ascii_strncasecmp(*buff,"miller",6) == 0)
    {
/* init space group (latmat? - if so, remove the make_latmat() in cell parse) */
    if (cflag && sflag)
      {
      if (!gflag)
        {
        if (space_lookup(data))
          printf("Error in space group lookup.\n");
        gflag = TRUE;
        }
      }
    else
      {
      if (data->id == MORPH)
        {
        printf("load_planes() error: miller encountered before space or cell.\n");
        return(2);
        }
      }

    if (num_tokens >= 4)
      {
      h = (gint) str_to_float(*(buff+1));
      k = (gint) str_to_float(*(buff+2));
      l = (gint) str_to_float(*(buff+3));
#if DEBUG_LOAD_PLANES
printf("read plane: %d %d %d\n", h, k, l);
#endif
      VEC3SET(m, h, k, l);
      plane = plane_find(m, data);
      if (!plane)
        {
        plane = plane_new(m, data);
        if (plane)
          new_planes = g_slist_prepend(new_planes, plane);
        }
      }
    }

/* potential data */
  if (num_tokens >= 8 && plane)
    {
/* NB: use create_shift(), as it sets some important defaults */
    shift = shift_new(0.0);
    if (shift)
      {
      shift->shift = str_to_float(*(buff+0));
      shift->region[0] = str_to_float(*(buff+1));
      shift->region[1] = str_to_float(*(buff+2));
      shift->esurf[0] = str_to_float(*(buff+3));
      shift->eatt[0] = str_to_float(*(buff+4));
      shift->esurf[1] = str_to_float(*(buff+5));
      shift->eatt[1] = str_to_float(*(buff+6));
      shift->gnorm = str_to_float(*(buff+7));
#if DEBUG_LOAD_PLANES
printf("adding shift: %f\n", shift->shift);
#endif
/* append to preserve order (eg import on existing plane set) */
      plane->shifts = g_slist_append(plane->shifts, shift);
      }
    }
  g_strfreev(buff);
  }
data->planes = g_slist_concat(data->planes, g_slist_reverse(new_planes));

/* compute dhkl's for the plane's list */
for (list=data->planes ; list ; list=g_slist_next(list))
  {
  plane = list->data;

/* create default shift if none found */
  if (!plane->shifts)
    {
    shift = shift_new(0.0);
    if (shift)
      plane->shifts = g_slist_append(plane->shifts, shift);
    }
/* get best energy for the plane */
  update_plane_energy(plane, data);
  }

/* compute symmetry related faces */
surf_symmetry_generate(data);

fclose(fp);
return(0);
}

/*********************************************/
/* determine if two hkl faces are equivalent */
/*********************************************/
#define DEBUG_FACET_EQUIV 0
gint facet_equiv(struct model_pak *data, gint *f1, gint *f2)
{
gint i, index[3];
gdouble vec[3];

/* NB: Dhkl test is just plain wrong, some symmetry related faces */
/* will have exactly the same Dhkl and be different surfaces */

#if DEBUG_FACET_EQUIV
printf("[%d %d %d] and (%d %d %d)", f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]);
printf(" [%.10f]\n", fabs(d2-d1));
#endif

/* symmetry test */
if (data->sginfo.spacenum > 0)
  {
  for (i=0 ; i<data->sginfo.order ; i++)
    {
    ARR3SET(vec, f1);
    vecmat(*(data->sginfo.matrix+i), vec);
    ARR3SET(index, vec);

#if DEBUG_FACET_EQUIV
P3VEC("op -> ", vec);
#endif

    if (index[0] == *f2 && index[1] == *(f2+1) && index[2] == *(f2+2))
      {
#if DEBUG_FACET_EQUIV
printf(" *** equivalent (symop %d).\n", i);
#endif
      return(TRUE);  
      }
    }
  }
else
  printf("No space group information.\n");

#if DEBUG_FACET_EQUIV
printf("not equivalent.\n");
#endif

return(FALSE);
}

/***********************************/
/* get a list of equivalent facets */
/***********************************/
#define DEBUG_GET_FACET_EQUIV 0
GSList *get_facet_equiv(struct model_pak *data, gint *f1)
{
gint i, flag, *f2, index[3];
gdouble d1, d2, vec[3];
GSList *flist=NULL, *list=NULL;

g_return_val_if_fail(data != NULL, NULL);
g_return_val_if_fail(f1 != NULL, NULL);

/* NEW */
ARR3SET(vec, f1);
vecmat(data->rlatmat, vec);
d1 = 1.0/VEC3MAG(vec);

#if DEBUG_GET_FACET_EQUIV
printf("search for equiv faces to (%d %d %d) (Dhkl = %f)\n",
        f1[0], f1[1], f1[2], d1);
#endif

/* add the supplied face to the list (needed to eliminate repetitions) */
f2 = g_malloc(3*sizeof(gint));
ARR3SET(f2, f1);
flist = g_slist_prepend(flist, (gpointer) f2);

if (data->sginfo.spacenum)
  {
/* skip 1st (trivial) op */
  for (i=1 ; i<data->sginfo.order ; i++)
    {
/* generate new symmetry related hkl */
    ARR3SET(vec, f1);
    vecmat(*(data->sginfo.matrix+i), vec);
    ARR3SET(index, vec);

/* NEW - weed out symop generated faces with different Dhkl values */
/* FIXME - why does this happen??? */
vecmat(data->rlatmat, vec);
d2 = 1.0/VEC3MAG(vec);

#if DEBUG_GET_FACET_EQUIV
printf("candidate: (%3d %3d %3d) : (Dhkl=%7.4f)", index[0], index[1], index[2], d2);
#endif

if (fabs(d2-d1) > FRACTION_TOLERANCE)
  {
#if DEBUG_GET_FACET_EQUIV
printf("[NO : dhkl mis-match]\n");
#endif
  continue;
  }

/* add new hkl if not found in the list */
    flag = 0;
    list = flist;
    while (list != NULL)
      {
      f2 = (gint *) list->data;
      if (index[0] == f2[0] && index[1] == f2[1] && index[2] == f2[2])
        {
        flag++;
        break;
        }
      list = g_slist_next(list);
      }
    if (!flag)
      {
      f2 = g_malloc(3*sizeof(gint));
      ARR3SET(f2, index);
      flist = g_slist_prepend(flist, f2);

#if DEBUG_GET_FACET_EQUIV
printf("[YES : symop %d]\n", i);
#endif
      }
#if DEBUG_GET_FACET_EQUIV
    else
      {
printf("[NO : already exists]\n");
      }
#endif
    }
  }
else
  printf("No space group information.\n");

flist = g_slist_reverse(flist);
return(flist);
}

/**************************************/
/* determine if systematically absent */
/**************************************/
#define DEBUG_FACET_ABSENT 0
gint surf_sysabs(struct model_pak *data, gint h, gint k, gint l)
{
gint i, flag;
gint f1[3] /*, f2[3]*/;
gdouble dp, dummy, vec[3], df[3];

#if DEBUG_FACET_ABSENT
printf("testing %d %d %d\n", h, k, l);
#endif

/* apply the symmetry matrices (skip identity) */
VEC3SET(f1, h, k, l);
for (i=1 ; i<data->sginfo.order ; i++)
  {
/* NEW - identity + translation alone is insufficent */
  if (matrix_is_identity(*(data->sginfo.matrix+i)))
    continue;

  VEC3SET(vec, h, k, l);

  vecmat(*(data->sginfo.matrix+i), vec);
  flag = 0;

/* not needed - we've expanded all symmetry operations (incl. inversion) */
/* test f1 . m = -f1 */
/*
  ARR3SET(df, f1);
  ARR3ADD(df, vec);
  if (VEC3MAGSQ(df) < FRACTION_TOLERANCE)
    if (data->sginfo.inversion)
      flag++;
*/
 
/* test f1 . m = f1 */
  ARR3SET(df, f1);
  ARR3SUB(df, vec);
  if (VEC3MAGSQ(df) < POSITION_TOLERANCE)
    flag++;

  if (flag)
    {
#if DEBUG_FACET_ABSENT
printf("symop [%d] satisfies 1st absence condition.\n", i);
P3MAT("matrix:", *(data->sginfo.matrix+i));
P3VEC("offset:", *(data->sginfo.offset+i));
#endif

/* test if <f1,t> = non-integer */
    ARR3SET(vec, *(data->sginfo.offset+i));
    ARR3MUL(vec, f1);
    dp = fabs(vec[0] + vec[1] + vec[2]);
    if (modf(dp, &dummy) > POSITION_TOLERANCE)
      {
#if DEBUG_FACET_ABSENT
printf("symop [%d] [%f %f %f] satisfies 2nd absence condition.\n",
       i, *(*(data->sginfo.offset+i)+0), *(*(data->sginfo.offset+i)+1),
          *(*(data->sginfo.offset+i)+2));
printf("facet is extinct.\n");
#endif
      return(TRUE);
      }
    }
  }

#if DEBUG_FACET_ABSENT
printf(">>> Not absent\n");
#endif

return(FALSE);
}

/***********************************/
/* determine if a plane is visible */
/***********************************/
#define DEBUG_FACET 0
gint facet_visible(struct model_pak *data, struct plane_pak *plane)
{
gdouble a, norm[3], vec[3];

/* NEW - an edge (2 vertices) is not to be treated as visible */
if (g_slist_length(plane->vertices) < 3)
  return(FALSE);

/* get plane  normal */
ARR3SET(norm, plane->norm);

VEC3MUL(norm, -1.0);

#if DEBUG_FACET
printf("facet (%2d%2d%2d), norm: %5.2f,%5.2f,%5.2f    " 
,plane->index[0] ,plane->index[1] ,plane->index[2]
,norm[0],norm[1],norm[2]);
#endif

/* absolute viewing vector, determining if visible or not */
camera_view(vec, data->camera);

/* got correct facet normal - is it visible? */
a = via(norm,vec,3);
if (a < 0.5*G_PI)
  {
#if DEBUG_FACET
printf("view: %5.2f,%5.2f,%5.2f (%5.1f) YES\n",vec[0],vec[1],vec[2],180.0*a/PI);
#endif
  return(TRUE);
  }
/* else... */
#if DEBUG_FACET
printf("view: %5.2f,%5.2f,%5.2f (%5.1f) NO\n",vec[0],vec[1],vec[2],180.0*a/PI);
#endif
return(FALSE);
}

/***********************/
/* Morphology creation */
/***********************/
#define DEBUG_READ_GMF 0
gint read_gmf(gchar *filename, struct model_pak *data)
{
/* checks */
g_assert(data != NULL);
g_assert(filename != NULL);

/* setup model parameters */
data->mode = FREE;
data->id = MORPH;
data->axes_type = OTHER;
data->morph_label = TRUE;
strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(data->filename);

/* read from a gmf file */
if (load_planes(filename, data))
  {
  printf("File load failed, check %s.\n", data->filename);
  model_delete(data);
  return(1);
  }

/* construct the hull (convert to H representation) */
data->num_vertices=0;
if (morph_build(data))
  return(3);

#if DEBUG_READ_GMF
printf("-------------------------\n");
printf("       Input planes: %d\n", data->num_planes);
printf("Calculated vertices: %d\n", data->num_vertices);
printf("-------------------------\n");
#endif

/* test success */
if (!g_slist_length(data->planes))
  return(4);

/* vertices are cartesian - so a fudge is needed */
model_prep(data);

/* proper rmax scaling & visibility calc */
coords_init(REDO_COORDS, data);

return(0);
}

