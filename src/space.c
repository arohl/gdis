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
#include <strings.h>
#include <ctype.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "matrix.h"
#include "sginfo.h"
#include "space.h"
#include "zone.h"
#include "dialog.h"
#include "interface.h"

/**********/
/* sginfo */
/**********/
static int str_ibegin(const char *s1, const char *s2) /* string ignore-case */
{                                                     /* begin              */
char u1, u2;

while (*s1 && *s2)
  {
  u1 = toupper(*s1++);
  u2 = toupper(*s2++);
  if (u1 < u2)
    return -1;
  else if (u1 > u2)
    return  1;
  }
if (*s2)
  return -1;
return 0;
}

/********************************/
/* cleanup raw sginfo structure */
/********************************/
void sginfo_free(T_SgInfo *SgInfo)
{
g_free(SgInfo->ListSeitzMx);
g_free(SgInfo->ListRotMxInfo);
}

/**********/
/* sginfo */
/**********/
/* TODO - we know what the input will be - trim */
gint BuildSgInfo(T_SgInfo *SgInfo, const gchar *SgName)
{
gint                VolLetter;
const T_TabSgName  *tsgn;


/* look for "VolA", "VolI", or "Hall" */

while (*SgName && isspace((int) *SgName)) SgName++;

VolLetter = -1;

if      (isdigit((int) *SgName))
  VolLetter = 'A';
else if (str_ibegin(SgName, "VolA") == 0)
  {
  VolLetter = 'A';
  SgName += 4;
  }
else if (str_ibegin(SgName, "VolI") == 0 || str_ibegin(SgName, "Vol1") == 0)
  {
  VolLetter = 'I';
  SgName += 4;
  }
else if (str_ibegin(SgName, "Hall") == 0)
  {
  VolLetter = 0;
  SgName += 4;
  }

while (*SgName && isspace((int) *SgName)) 
  SgName++;

/* default is "VolA" */

if (VolLetter == -1)
  VolLetter = 'A';

/* if we don't have a Hall symbol do a table look-up */

tsgn = NULL;

if (VolLetter)
  {
  tsgn = FindTabSgNameEntry(SgName, VolLetter);
  if (tsgn == NULL) return -1; /* no matching table entry */
  SgName = tsgn->HallSymbol;
  }

/* Allocate memory for the list of Seitz matrices and
   a supporting list which holds the characteristics of
   the rotation parts of the Seitz matrices
*/

SgInfo->MaxList = 192; /* absolute maximum number of symops */

SgInfo->ListSeitzMx
    = g_malloc(SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));

SgInfo->ListRotMxInfo
    = g_malloc(SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));

/* Initialize the SgInfo structure */

InitSgInfo(SgInfo);
SgInfo->TabSgName = tsgn; /* in case we know the table entry */

/* Translate the Hall symbol and generate the whole group */

ParseHallSymbol(SgName, SgInfo);
if (SgError != NULL)
  return -1;

  /* Do some book-keeping and derive crystal system, point group,
     and - if not already set - find the entry in the internal
     table of space group symbols
   */

return CompleteSgInfo(SgInfo);
}

/*****************************************/
/* generate positions of a complete cell */
/*****************************************/
void space_fill_cell(struct model_pak *model)
{
gint j, s, smin;
GSList *cores=NULL, *list;
struct core_pak *core, *sr_core;
struct shel_pak *sr_shell;

/* include an inversion (ie mult by -1) operation? */
smin = 0;
if (model->sginfo.inversion == -1)
  smin = -2;

/* s = +1 or -1 for inversion of matrix elements */
for (s=1 ; s>smin ; s-=2)
  {
/* loop over group operations */
/* NB: include j=0 in case of inversion centre ie (-x,-y,-z) */
  for (j=(s+1)/2 ; j<model->sginfo.order ; j++)
    {
/* loop over asymetric cores */
    for (list=model->cores ; list ; list=g_slist_next(list))
      {
      core = list->data;
      if (!core->primary)
        continue;

/* add a symmetry related core */
      sr_core = dup_core(core);
      cores = g_slist_prepend(cores, sr_core);
      sr_core->primary = FALSE;
      sr_core->primary_core = core;
      VEC3MUL(sr_core->x, s);
      vecmat(*(model->sginfo.matrix+j), sr_core->x);
      sr_core->x[0] += *(*(model->sginfo.offset+j)+0);
      sr_core->x[1] += *(*(model->sginfo.offset+j)+1);
      sr_core->x[2] += *(*(model->sginfo.offset+j)+2);

/* add a symmetry related shell */
      if (sr_core->shell)
        {
        sr_shell = sr_core->shell;
        model->shels = g_slist_prepend(model->shels, sr_shell);
        sr_shell->primary = FALSE;
        sr_shell->primary_shell = sr_shell;
        VEC3MUL(sr_shell->x, s);
        vecmat(*(model->sginfo.matrix+j), sr_shell->x);
        sr_shell->x[0] += *(*(model->sginfo.offset+j)+0);
        sr_shell->x[1] += *(*(model->sginfo.offset+j)+1);
        sr_shell->x[2] += *(*(model->sginfo.offset+j)+2);
        }
      }
    }
  }

/* update lists */
model->cores = g_slist_concat(model->cores, g_slist_reverse(cores));

/* remove any duplicates */
/* NB: do this before pbc constrain, as it updates core-shell linkages */
zone_init(model);
delete_duplicate_cores(model);
shell_make_links(model);
}

/****************************************/
/* initialize the space group structure */
/****************************************/
void space_init(gpointer data)
{
struct space_pak *space=data;

g_assert(space != NULL);

space->lookup=TRUE;
space->spacenum=0;
space->lattice=0;
space->pointgroup=0;
space->cellchoice=0;
space->inversion=FALSE;
space->order=0;
space->spacename=NULL;
space->latticename=NULL;
space->matrix=NULL;
space->offset=NULL;
}

/********************************/
/* free a space group structure */
/********************************/
void space_free(gpointer data)
{
gint i;
struct space_pak *space=data;

g_assert(space != NULL);

/* free data */
g_free(space->spacename);
g_free(space->latticename);
for (i=space->order ; i-- ; )
  {
  g_free(*(space->matrix+i));
  g_free(*(space->offset+i));
  }
g_free(space->matrix);
g_free(space->offset);

/* enforce blankness */
space_init(space);
}

/********************************/
/* hexagonal/rhombohedral check */
/********************************/
gint space_primitive_cell(struct model_pak *model)
{
gdouble dx, x2;

dx = (model->pbc[0] - model->pbc[1]);
x2 = dx*dx;
dx = (model->pbc[0] - model->pbc[2]);
x2 += dx*dx;
dx = (model->pbc[1] - model->pbc[2]);
x2 += dx*dx;

/* better tolerances? */
if (x2 < 0.001)
  {
/*
printf("Treating as rhombohedral.\n");
*/
  return(TRUE);
  }
/*
printf("Treating as hexagonal.\n");
*/
return(FALSE);
}

/**************************/
/* do space group loopkup */
/**************************/
#define DEBUG_GEN_POS 0
gint space_lookup(struct model_pak *data)
{
gint f, i, j, nt;
gint m;
gint iList, nTrV, iTrV, nLoopInv, iLoopInv, iMatrix;
gchar *label;
GString *name;
const T_RTMx  *lsmx;
const gint *r, *t, *TrV;
T_RTMx SMx;
T_SgInfo SgInfo;

g_return_val_if_fail(data != NULL, 1);

/* NEW - can flag to skip this eg GULP animation */
/* where symmetry is broken in subsequent frames */
if (!data->sginfo.lookup)
  {
  return(0);
  }

/* request info via number or name */
/* cell choice (TODO - neaten this ugly code) */
name = g_string_new(NULL);
if (data->sginfo.spacename)
  g_string_sprintf(name, "%s", g_strstrip(data->sginfo.spacename));
else
  {
  if (data->sginfo.spacenum > 0)
    g_string_sprintf(name, "%d", data->sginfo.spacenum);
  else
    g_string_sprintf(name, "P 1");
  }

/* cell choice */
if (data->sginfo.cellchoice)
  g_string_sprintfa(name, ":%d", data->sginfo.cellchoice);

/* if trigonal lattice, determine if cell is in hex or rhombo format */
if (g_strrstr(name->str, "R") || data->sginfo.spacenum == 146
                              || data->sginfo.spacenum == 148
                              || data->sginfo.spacenum == 155
                              || data->sginfo.spacenum == 160
                              || data->sginfo.spacenum == 161
                              || data->sginfo.spacenum == 166
                              || data->sginfo.spacenum == 167)
  {
  if (space_primitive_cell(data))
    g_string_sprintfa(name, ":R");
  }

#if DEBUG_GEN_POS
printf("retrieving as [%s]\n", name->str);
#endif

/* main call */
if (BuildSgInfo(&SgInfo, name->str) != 0)
  {
/* CURRENT - if order > 0 -> fill the cell */
  if (data->sginfo.order)
    {
    gui_text_show(WARNING, "No spacegroup: using explicit matrices.\n");
    space_fill_cell(data);
    return(0);
    }
  return(2);
  }
g_string_free(name, TRUE);

/* process returned space group name */
label = g_strdup(SgInfo.TabSgName->SgLabels);
for (i=0 ; i<strlen(label) ; i++)
  if (*(label+i) == '_')
    *(label+i) = ' ';

/* copy label to avoid any funny GULP characters */
i=strlen(label);
while(i>0)
  {
  if (*(label+i) == '=')
    {
    i++;
    break;
    }
  i--;
  }

/* fill in space group name */
if (!data->sginfo.spacename)
  data->sginfo.spacename = g_strdup_printf("%s", label+i);
g_free(label);

/* fill in number (if an invalid value was supplied) */
if (data->sginfo.spacenum < 1)
  data->sginfo.spacenum = SgInfo.TabSgName->SgNumber;

/* crystal lattice type */
data->sginfo.lattice = SgInfo.XtalSystem;

if (data->sginfo.latticename)
  g_free(data->sginfo.latticename);

switch (SgInfo.XtalSystem)
  {
  case XS_CUBIC:
    data->sginfo.latticename = g_strdup("Cubic");
    break;
  case XS_TETRAGONAL:
    data->sginfo.latticename = g_strdup("Tetragonal");
    break;
  case XS_ORTHOROMBIC:
    data->sginfo.latticename = g_strdup("Orthorhombic");
    break;
  case XS_HEXAGONAL:
    data->sginfo.latticename = g_strdup("Hexagonal");
    break;
  case XS_TRIGONAL:
    data->sginfo.latticename = g_strdup("Trigonal");
    break;
  case XS_MONOCLINIC:
    data->sginfo.latticename = g_strdup("Monoclinic");
    break;
  case XS_TRICLINIC:
    data->sginfo.latticename = g_strdup("Triclinic");
    break;
  default:
    data->sginfo.latticename = g_strdup("Unknown");
  }

/* point group */
data->sginfo.pointgroup = SgInfo.PointGroup;
data->sginfo.inversion = SgInfo.Centric;
data->sginfo.centering = SgInfo.LatticeInfo->Code;

#if DEBUG_GEN_POS
printf("   Space group name: %s\n", data->sginfo.spacename);
printf("        Point group: %d\n", data->sginfo.pointgroup);
printf(" Space group number: %d\n", data->sginfo.spacenum);
printf("       Lattice type: %c\n", data->sginfo.centering);
printf("              Order: %d\n", SgInfo.OrderL);
printf("   Inversion center: ");
if (data->sginfo.inversion == -1)
  printf("yes\n");
else
  printf("no\n");
#endif

/* CURRENT */
if (!data->sginfo.order)
  {
/* allocate for order number of pointers (to matrices) */
data->sginfo.order = SgInfo.OrderL;
data->sginfo.matrix = (gdouble **) g_malloc(SgInfo.OrderL*sizeof(gdouble *));
data->sginfo.offset = (gdouble **) g_malloc(SgInfo.OrderL*sizeof(gdouble *));

iMatrix = 0;

nLoopInv = Sg_nLoopInv(&SgInfo);

nTrV = SgInfo.LatticeInfo->nTrVector;
 TrV = SgInfo.LatticeInfo->TrVector;

/* matrix counter */
m=0;
for (iTrV = 0; iTrV < nTrV; iTrV++, TrV += 3)
  {
  for (iLoopInv = 0; iLoopInv < nLoopInv; iLoopInv++)
    {
    if (iLoopInv == 0)
      f =  1;
    else
      f = -1;

    lsmx = SgInfo.ListSeitzMx;

/* loop over all matrices (order of the group) */
    for (iList = 0; iList < SgInfo.nList; iList++, lsmx++)
      {
      for (i = 0; i < 9; i++)
        SMx.s.R[i] = f * lsmx->s.R[i];
      for (i = 0; i < 3; i++)
        SMx.s.T[i] = iModPositive(f * lsmx->s.T[i] + TrV[i], STBF);

      r = SMx.s.R;
      t = SMx.s.T;

      *(data->sginfo.matrix+m) = (gdouble *) g_malloc(9*sizeof(gdouble));
      *(data->sginfo.offset+m) = (gdouble *) g_malloc(3*sizeof(gdouble));

      for (i = 0; i < 3; i++, t++)
        {
        for (j = 0; j < 3; j++, r++)
          {
/* mth general position */
          *(*(data->sginfo.matrix+m)+3*i+j) = (gdouble) *r;
          }
        nt = iModPositive(*t, STBF);
        if (nt >  STBF / 2)
          nt -= STBF;
/* offset matrix (3x1) */
        *(*(data->sginfo.offset+m)+i) = (gdouble) nt / (gdouble) STBF;
        }
      m++;
      }
    }
  }

/* number of matrices matches the order? */
  if (m != data->sginfo.order)
    printf("Serious error in space_lookup()\n");
  }
else
  {
#if DEBUG_GEN_POS
printf("Skipping symmetry matrix generation.\n");
#endif
  m = data->sginfo.order;
  }

#if DEBUG_GEN_POS
printf("Symmetry matrices: %d\n", m);
for (i=0 ; i<m ; i++)
  {
  printf("matrix %d",i);

P3MAT(" : ", *(data->sginfo.matrix+i));

/*
  for (j=0 ; j<9 ; j++)
    printf(" %4.1f",*(*(data->sginfo.matrix+i)+j));
  printf("\n");
*/

  printf("offset %d :",i);
  for (j=0 ; j<3 ; j++)
    printf(" %4.1f",*(*(data->sginfo.offset+i)+j));
  printf("\n");
  }
#endif

/* deprec - raw sginfo structure */
/*
if (data->sginfo.raw)
  g_free(data->sginfo.raw);
data->sginfo.raw = g_malloc(sizeof(SgInfo));
memcpy(data->sginfo.raw, &SgInfo, sizeof(SgInfo));
*/
sginfo_free(&SgInfo);

/* generate symmetry related cores */
space_fill_cell(data);

#if DEBUG_GEN_POS
printf("[*] num atoms: %d\n", g_slist_length(data->cores));
printf("[*] num shels: %d\n", g_slist_length(data->shels));

#endif

return(0);
}

/************************************/
/* remove all symmetry from a model */
/************************************/
void space_make_p1(struct model_pak *model)
{
GSList *list;
struct core_pak *core;
struct shel_pak *shell;

/* checks */
g_assert(model != NULL);

/* initialize atoms */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->orig = TRUE;
  core->primary = TRUE;
  }
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shell = list->data;
  shell->orig = TRUE;
  shell->primary = TRUE;
  }

/* set spacegroup */
space_free(&model->sginfo);
model->sginfo.spacenum = 1;
model->sginfo.cellchoice = 0;
space_lookup(model);
}

/***********************************************/
/* convert periodic images into full supercell */
/***********************************************/
void space_make_supercell(struct model_pak *model)
{
gint i;
gdouble v[3], m[3];
GSList *list, *ilist, *clist, *slist;
struct core_pak *core, *copy;
struct shel_pak *shell;
struct image_pak *image;

/* checks */
if (!model)
  return;

/* all (non-deleted) atoms are now primary & non fractional */
delete_commit(model);
clist = slist = NULL;

/* loop over images 1st, so we preserve core ordering in cell images */
for (ilist=model->images ; ilist ; ilist=g_slist_next(ilist))
  {
  image = ilist->data;
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    copy = dup_core(core);
    ARR3ADD(copy->x, image->pic);
    clist = g_slist_prepend(clist, copy);
    if (copy->shell)
      {
      shell = copy->shell;
      ARR3ADD(shell->x, image->pic);
      slist = g_slist_prepend(slist, shell);
      }
    }
  }
free_slist(model->images);
model->images = NULL;

/* order preserving core/shell additions */
model->cores = g_slist_concat(model->cores, g_slist_reverse(clist));
model->shels = g_slist_concat(model->shels, g_slist_reverse(slist));

/* get multiple of lattice vectors desired */
m[0] = fabs(model->image_limit[1] + model->image_limit[0]);
m[1] = fabs(model->image_limit[3] + model->image_limit[2]);
m[2] = fabs(model->image_limit[5] + model->image_limit[4]);

/* scale up the energies */
model->gulp.energy *= m[0] * m[1] * m[2];
model->gulp.sbulkenergy *= m[0] * m[1] * m[2];

/* scale the fractional coordinates down */
model->fractional = TRUE;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  core->x[0] /= m[0];
  core->x[1] /= m[1];
  core->x[2] /= m[2];
  }
for (list=model->shels ; list ; list=g_slist_next(list))
  {
  shell = list->data;
  shell->x[0] /= m[0];
  shell->x[1] /= m[1];
  shell->x[2] /= m[2];
  }

/* NB: need to scale the pbc so zone init works */
ARR3MUL(model->pbc, m);

/* get current lattice vectors and scale up */
for (i=0 ; i<3 ; i++)
  {
  v[0] = model->latmat[i];
  v[1] = model->latmat[i+3];
  v[2] = model->latmat[i+6];
  VEC3MUL(v, m[i]);
  model->latmat[i] = v[0];
  model->latmat[i+3] = v[1];
  model->latmat[i+6] = v[2];
  }

/* update cell param related data */
model->construct_pbc = TRUE;

/* TODO - currently, tell gdis to completely redo pbond calc again */
/* TODO - a more efficient way would be to recalc */
model->done_pbonds = FALSE;

model->image_limit[0] = 0.0;
model->image_limit[1] = 1.0;
model->image_limit[2] = 0.0;
model->image_limit[3] = 1.0;
model->image_limit[4] = 0.0;
model->image_limit[5] = 1.0;

space_make_p1(model);
}

/***********************************/
/* periodic image creation routine */
/***********************************/
#define DEBUG_UPDATE_IMAGES 0
void space_make_images(gint mode, struct model_pak *model)
{
gint a, b, c, i;
gint num_cells, num_images;
struct image_pak *image;

/* checks */
g_assert(model != NULL);

#if DEBUG_UPDATE_IMAGES
printf("---------------------------\n");
printf("model : %d\n", model->number);
printf("a : %f - %f\n", model->image_limit[0], model->image_limit[1]);
printf("b : %f - %f\n", model->image_limit[2], model->image_limit[3]);
printf("c : %f - %f\n", model->image_limit[4], model->image_limit[5]);
printf("---------------------------\n");
#endif

num_images = 0;
/* check mode */
switch(mode)
  {
  case CREATE:
/* update image list */
    free_slist(model->images);
    model->images = NULL;
    num_cells = (model->image_limit[0]+model->image_limit[1])
              * (model->image_limit[2]+model->image_limit[3])
              * (model->image_limit[4]+model->image_limit[5]);
    num_images = num_cells-1;
/* setup for pic iteration */
    a = -model->image_limit[0];
    b = -model->image_limit[2];
    c = -model->image_limit[4];
    for (i=num_cells ; i-- ; )
      {
/* image increment */
      if (a == model->image_limit[1])
        {
        a = -model->image_limit[0];
        b++;
        if (b == model->image_limit[3])
          {
          b = -model->image_limit[2];
          c++;
          if (c == model->image_limit[5])
            break;
          }
        }
/* don't duplicate the primary cell */
      if (a || b || c)
        {
        image = g_malloc(sizeof(struct image_pak));

        VEC3SET(image->pic, a, b, c);
        model->images = g_slist_prepend(model->images, image);
        }
      a++;
      }
    break;

  case INITIAL:
/* reset values to default */
    for (i=0 ; i<model->periodic ; i++)
      {
      model->image_limit[2*i] = 0;
      model->image_limit[2*i+1] = 1;
      }
/* delete all images */
    free_slist(model->images);
    model->images = NULL;
    break;
  }

/*
gui_relation_update(model);
*/
}

