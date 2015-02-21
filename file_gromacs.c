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
#include "matrix.h"
#include "measure.h"
#include "model.h"
#include "parse.h"
#include "file.h"
#include "render.h"
#include "scan.h"
#include "select.h"
#include "ff.h"
#include "ff_pak.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/************************************************/
/* convert internal forcefield type to GROMACS  */
/************************************************/
gpointer gromacs_type(struct forcefield_pak *ff)
{
struct forcefield_pak *gff;

gff = ff_dup(ff);

switch (ff->type)
  {
/* bonds */
  case FF_HARMONIC:
    switch (ff->data_expected)
       {
/* harmonic */
       case 1:
         gff->type = 1;
         gff->data_expected = 1;
         break;
/* cubic */
       case 2:
       case 3:
         gff->type = 4;
         gff->data_expected = 2;
         break;

       default:
         g_assert_not_reached();
       }
    break;

/* morse */
  case FF_MORSE:
    gff->type = 2;
    gff->data_expected = 2;
    break;

/* three body */
  case FF_3B_HARMONIC:
    switch (ff->data_expected)
      {
      case 3:
        gff->type = 3;
        gff->data_expected = 5;
        break;

      case 1:
      case 2:
      default:
        gff->type = 1;
        gff->data_expected = 1;
       }
    break;

/* four body */
  case FF_DIHEDRAL:
    gff->type = 1;
    gff->data_expected = 1;
    break;
  case FF_IMPROPER:
    gff->type = 2;
    gff->data_expected = 1;
    break;
  case FF_DIHEDRAL_RB:
    gff->type = 3;
    gff->data_expected = 6;
    break;

/* non bond */
  case FF_LENNARD:
    gff->type = 1;
    gff->data_expected = 2;
    break;
  case FF_BUCKINGHAM:
    gff->type = 2;
    gff->data_expected = 2;
    break;
  }

/* GROMACS length units are nm */
switch (ff->bond_units)
  {
  case FF_ANG:
    gff->bond_value *= 0.1;
    break;
  }

return(gff);
}

/******************************/
/* raw coordinate data (.gro) */
/******************************/
gint read_gromacs_gro(gchar *filename, struct model_pak *model)
{
gint num_tokens;
gchar *line, **buff;
FILE *fp;

g_assert(model != NULL);

fp = fopen(filename,"rt");
if (!fp)
  return(1);

/* skip title */
line = file_read_line(fp);
g_free(line);
line = file_read_line(fp);

while (line)
  {

/* FIXME - properly should do formatted read, since things can be */
/* adjacent with no spaces in between -> tokenize() will fail */
/*
  fprintf(fp, "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
               n, "UNK", core->atom_label, 1, 
               x[0], x[1], x[2], 0.1*core->v[0], 0.1*core->v[1], 0.1*core->v[2]);
*/

  if (strlen(line) > 10)
    {
/* tokenize from atom label (avoids some problems with no space separation) */
    buff = tokenize(&line[10], &num_tokens);

    if (num_tokens > 4)
      {
/* TODO - deal with water properly (ie set atom_types) */
      if (g_ascii_strncasecmp(*(buff+0), "HW", 2) == 0)
        {
        g_free(*(buff+0));
        *(buff+0) = g_strdup("H");
        }
      if (g_ascii_strncasecmp(*(buff+0), "OW", 2) == 0)
        {
        g_free(*(buff+0));
        *(buff+0) = g_strdup("O");
        }

      if (elem_test(*(buff+0)))
        {
        struct core_pak *core = new_core(*(buff+0), model);

        core->x[0] = 10.0*str_to_float(*(buff+2));
        core->x[1] = 10.0*str_to_float(*(buff+3));
        core->x[2] = 10.0*str_to_float(*(buff+4));

        model->cores = g_slist_prepend(model->cores, core);
        }
      }
    g_strfreev(buff);
    }
  g_free(line);
  line = file_read_line(fp);
  }

model->cores = g_slist_reverse(model->cores);

model_prep(model);

return(0);
}

/******************************/
/* raw coordinate data (.gro) */
/******************************/
gint write_gromacs_gro(gchar *filename, struct model_pak *model)
{
gint n, c, i, s;
gdouble x[3], min[3], max[3];
GSList *list;
struct core_pak *core;
FILE *fp;

fp = fopen(filename, "wt");
if (!fp)
  return(1);

/* init min/max for box calculation */
for (i=3 ; i-- ; )
  {
  max[i] = -G_MAXDOUBLE;
  min[i] = G_MAXDOUBLE;
  }

fprintf(fp, "GROMACS coordinates (generated by GDIS v%f)\n", VERSION);

c = g_slist_length(model->cores);
s = g_slist_length(model->shels);

fprintf(fp, "%5d\n", c+s);

n = 1;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* cartesians */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
/* angs -> nm */
  VEC3MUL(x, 0.1);
  fprintf(fp, "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
               n, "UNK", core->atom_label, 1, 
               x[0], x[1], x[2], 0.1*core->v[0], 0.1*core->v[1], 0.1*core->v[2]);

/* record min/max for box calculation */
  for (i=3 ; i-- ; )
    {
    if (x[i] < min[i])
      min[i] = x[i];
    if (x[i] > max[i])
      max[i] = x[i];
    }

  if (core->shell)
    {
/* TODO - has an 'x' post-fix to the atom_label */
    }

  n++;
  }

/* this is the distance between the min and max atom coords in x, y, and z */
/* box calculation */
if (model->periodic == 3)
  {
  fprintf(fp, "% f %f %f\n", 0.1*model->pbc[0], 0.1*model->pbc[1], 0.1*model->pbc[2]);
  }
else
  {
  for (i=0 ; i<3 ; i++)
    fprintf(fp, " %f", fabs(max[i] - min[i]));
  fprintf(fp, "\n");
  }

fclose(fp);

return(0);
}

/********************/
/* write the header */
/********************/
gint write_gromacs_header(FILE *fp, struct model_pak *model)
{

fprintf(fp, "\n;\n; Generated by GDIS v%f\n;\n", VERSION);

/* put in footer? */
fprintf(fp, "\n[ system ]\n");
fprintf(fp, "%s\n", model->basename);

fprintf(fp, "\n[ moleculetype ]\n");
/* label, exclusion type */
fprintf(fp, "%s %d\n", model->basename, 3);

return(0);
}

/********************/
/* write the footer */
/********************/
gint write_gromacs_footer(FILE *fp, struct model_pak *model)
{

fprintf(fp, "\n[ molecules ]\n");
/* label, number of occurences */
fprintf(fp, "%s 1\n", model->basename);

return(0);
}

/*******************/
/* write the atoms */
/*******************/
gint write_gromacs_atoms(FILE *fp, struct model_pak *model)
{
gint n;
gdouble w;
GSList *list;
struct core_pak *core;

fprintf(fp, "\n[ atoms ]\n");

n = 1;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  w = elements[core->atom_code].weight;

/* if 2nd column is atom type - use if available */
  if (core->atom_type)
    {
    fprintf(fp, "%9d %5s %5d %5s %5s %5d %8.4f %8.4f",
                 n, core->atom_type, 1, "UNK", core->atom_label, n,
                 core->charge, w);
    }
  else
    {
    fprintf(fp, "%9d %5s %5d %5s %5s %5d %8.4f %8.4f",
                 n, core->atom_label, 1, "UNK", core->atom_label, n,
                 core->charge, w);
    }

/* the lambda = 1.0 values - default to dummy for the time being */
  fprintf(fp, "    DUM  0.0  %8.4f\n", w);

  n++;
  }

return(0);
}

/********************************/
/* write forcefield information */
/********************************/
gint write_gromacs_bonds(FILE *fp, struct model_pak *model)
{
gint i, j, n;
GSList *list;
struct core_pak *core[2];
struct forcefield_pak *ff, *ff_gromacs;

g_assert(model != NULL);
g_assert(fp != NULL);

/* header */
fprintf(fp, "\n[ bonds ]\n");

/* for each bond FF type in the list -> ff_bond_search -> output */
for (list=model->bonds ; list ; list=g_slist_next(list))
  {
  struct bond_pak *bond = list->data;

  core[0] = bond->atom1;
  core[1] = bond->atom2;

  ff = ff_search(core, 2, model->ff_list);

  if (ff)
    {
/* get core numbering order */
    i = g_slist_index(model->cores, core[0]) + 1;
    j = g_slist_index(model->cores, core[1]) + 1;

    ff_gromacs = gromacs_type(ff);

/* lambda 0 */
/* print atoms and type */
    fprintf(fp, "%2d %2d   %d  ", i, j, ff_gromacs->type);

/* NB: gromacs has the bond length 1st, followed by the potential parameters */
/* which is the reverse from internal/GULP format */
    if (ff_gromacs->bond_expected)
      fprintf(fp, "%f", ff_gromacs->bond_value);
    for (n=0 ; n<ff_gromacs->data_expected ; n++)
      fprintf(fp, " %f", ff_gromacs->data[n]);

/* lambda 1 */
    fprintf(fp, "    ");
    if (ff_gromacs->bond_expected)
      fprintf(fp, "%f", ff_gromacs->bond_value);
    for (n=0 ; n<ff_gromacs->data_expected ; n++)
      fprintf(fp, " %f", ff_gromacs->data[n]);

    fprintf(fp, "\n");

    g_free(ff_gromacs);
    }
  } 

return(0);
}

/********************************/
/* write forcefield information */
/********************************/
#define MAX_CORES 10
gint write_gromacs_angles(FILE *fp, struct model_pak *model)
{
gint i, j, k, m, n, ni, nj;
GSList *list1, *list2, *nlist;
struct core_pak *core[3], *nc[MAX_CORES];
struct forcefield_pak *ff, *ff_gromacs;

g_assert(model != NULL);
g_assert(fp != NULL);

/* header */
fprintf(fp, "\n[ angles ]\n");

/* TODO - for all atoms with 2 + bonds -> check if 3 body term that matches */
for (list1=model->cores,k=1 ; list1 ; list1=g_slist_next(list1),k++)
  {
  core[1] = list1->data;

  nlist = connect_neighbours(core[1]);

  n = g_slist_length(nlist);

/* skip if too few, or suspiciously many */
  if (n < 2 || n >= MAX_CORES)
    continue;

  i = 0;
  for (list2 = nlist ; list2 ; list2=g_slist_next(list2))
    nc[i++] = list2->data;

/* iterate over unique neighbour pairs */
  for (i=0 ; i<n-1 ; i++)
    {
    core[0] = nc[i];

    for (j=i+1 ; j<n ; j++)
      {
      core[2] = nc[j];

      ff = ff_search(core, 3, model->ff_list);

      if (ff)
        {
/* TODO - write the entry */
/*
printf("FF found\n");
*/

/* get core numbering order */
        ni = g_slist_index(model->cores, core[0]) + 1;
        nj = g_slist_index(model->cores, core[2]) + 1;

        ff_gromacs = gromacs_type(ff);

/* lambda 0 */
/* print atoms and type */
        fprintf(fp, "%2d %2d %2d   %d  ", ni, k, nj, ff_gromacs->type);

/* NB: gromacs has the bond length 1st, followed by the potential parameters */
/* which is the reverse from internal/GULP format */
        if (ff_gromacs->bond_expected)
          fprintf(fp, "%f", ff_gromacs->bond_value);
        for (m=0 ; m<ff_gromacs->data_expected ; m++)
          fprintf(fp, " %f", ff_gromacs->data[m]);

/* lambda 1 */
        fprintf(fp, "    ");
        if (ff_gromacs->bond_expected)
          fprintf(fp, "%f", ff_gromacs->bond_value);
        for (m=0 ; m<ff_gromacs->data_expected ; m++)
          fprintf(fp, " %f", ff_gromacs->data[m]);

        fprintf(fp, "\n");
        }
      }
    }
  }

return(0);
}

/********************************/
/* write forcefield information */
/********************************/
gint write_gromacs_dihedrals(FILE *fp, struct model_pak *model)
{
gint i, j, k, l, n;
gdouble angle, da;
GSList *list, *list1, *list2, *end1, *end2;
GSList *d_list, *r_list, *i_list;
struct core_pak *c[4];
struct bond_pak *bond;
struct forcefield_pak *ff, *ff_gromacs;

g_assert(model != NULL);

d_list = ff_filter_list(FF_DIHEDRAL, model->ff_list);
r_list = ff_filter_list(FF_DIHEDRAL_RB, model->ff_list);
i_list = ff_filter_list(FF_IMPROPER, model->ff_list);

/* header */
fprintf(fp, "\n[ dihedrals ]\n");

/* generate torsions for all combinations about a bond */
for (list=model->bonds ; list ; list=g_slist_next(list))
  {
  bond = list->data;

  if (bond->type == BOND_HBOND)
    continue;

  c[1] = bond->atom1;
  c[2] = bond->atom2;

/* get lists of connected atoms at each end */
  end1 = connect_neighbours(c[1]);
  end1 = g_slist_remove(end1, c[2]);
  end2 = connect_neighbours(c[2]);
  end2 = g_slist_remove(end2, c[1]);

/* skip bonds with no neighbour atoms */
/*
  if (!end1 || !end2)
    goto gromacs_dihedral_loop_done;
*/

/* enumerate four body terms */
  for (list1=end1 ; list1 ; list1=g_slist_next(list1))
    {
    c[0] = list1->data;

    for (list2=end2 ; list2 ; list2=g_slist_next(list2))
      {
      c[3] = list2->data;

/* search for standard dihedral potential */
      ff = ff_search(c, 4, d_list);

/* search for RB dihedral */
      if (!ff)
        ff = ff_search(c, 4, r_list);

/* CURRENT */
/*
printf("[%s][%s][%s][%s]", c[0]->atom_label, c[1]->atom_label, c[2]->atom_label, c[3]->atom_label);
printf(" : dihedral = %f\n", measure_torsion(c));
*/

      if (ff)
        {
        ff_gromacs = gromacs_type(ff);

        i = g_slist_index(model->cores, c[0])+1;
        j = g_slist_index(model->cores, c[1])+1;
        k = g_slist_index(model->cores, c[2])+1;
        l = g_slist_index(model->cores, c[3])+1;

/* lambda 0 */
/* print atoms and type */
        fprintf(fp, "%2d %2d %2d %2d   %d  ", i, j, k, l, ff_gromacs->type);
/* NB: gromacs has the bond length 1st, followed by the potential parameters */
/* which is the reverse from internal/GULP format */
        if (ff_gromacs->bond_expected)
          fprintf(fp, "%f", ff_gromacs->bond_value);
        for (n=0 ; n<ff_gromacs->data_expected ; n++)
          fprintf(fp, " %f", ff_gromacs->data[n]);

/* lambda 1 */
        fprintf(fp, "    ");
        if (ff_gromacs->bond_expected)
          fprintf(fp, "%f", ff_gromacs->bond_value);
        for (n=0 ; n<ff_gromacs->data_expected ; n++)
          fprintf(fp, " %f", ff_gromacs->data[n]);

        fprintf(fp, "\n");

        g_free(ff_gromacs);
        }
      }
    }

  g_slist_free(end1);
  g_slist_free(end2);
  }

/* CURRENT - improper torsions (implemented as proper) */
/* FIXME - only valid for OPLS-AA */
/*
fprintf(fp, "\n[ impropers ]\n");
*/
fprintf(fp, "; Improper terms, implemented as proper torsions.\n");

/* for each atom with 3 bonds (planar) */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  c[0] = list->data;

  list1 = connect_neighbours(c[0]); 

/* CHECK - only an atom with 3 bonds can be a candidate */
  if (g_slist_length(list1) == 3)
    {
    c[1] = list1->data;
    list2 = g_slist_next(list1);
    c[2] = list2->data;
    list2 = g_slist_next(list2);
    c[3] = list2->data;


/* TODO - want c[3] to be defined so that the c[1] - c[2] axis is */
/* orthogonal to the c[0] - c[3] direction */

/* TODO - choose c[1] & c[2] to be same type (if possible) */

/* improper torsion - only allow if angle is ~180 */
    angle = measure_torsion(c);
/*
    da = fabs(180.0 - fabs(angle));
*/
    da = fabs(angle);

/* 10 degree tolerance */
    ff = NULL;
    if (da < 10.0)
      ff = ff_search(c, 4, i_list);

    if (ff)
      {
      ff_gromacs = gromacs_type(ff);

      i = g_slist_index(model->cores, c[0])+1;
      j = g_slist_index(model->cores, c[1])+1;
      k = g_slist_index(model->cores, c[2])+1;
      l = g_slist_index(model->cores, c[3])+1;

/* lambda 0 */
      fprintf(fp, "%2d %2d %2d %2d   1 ", i, j, k, l);
      fprintf(fp, " %f", ff_gromacs->bond_value);
      for (n=0 ; n<ff_gromacs->data_current ; n++)
        fprintf(fp, " %f", ff_gromacs->data[n]);

/* lambda 1 */
      fprintf(fp, "     %f", ff_gromacs->bond_value);
      for (n=0 ; n<ff_gromacs->data_current ; n++)
        fprintf(fp, " %f", ff_gromacs->data[n]);

      fprintf(fp, "\n");

      g_free(ff_gromacs);
      }
    }
  g_slist_free(list1);
  }

/* cleanup */
g_slist_free(d_list);
g_slist_free(r_list);
g_slist_free(i_list);

return(0);
}

/********************************/
/* write forcefield information */
/********************************/
gint write_gromacs_nonbond(FILE *fp, struct model_pak *model)
{
gint i;
GSList *list;
struct forcefield_pak *ff, *ff_gromacs;

g_assert(model != NULL);

/* header */
fprintf(fp, "\n[ nonbond_params ]\n");

for (list=model->ff_list ; list ; list=g_slist_next(list))
  {
  ff = list->data;

/* skip anything that doesnt have only 2 atom identifiers */
  if (ff->atoms_expected != 2)
    continue;
  if (ff->atoms_current != 2)
    continue;

/* GROMACS only supports LJ/BUCK */
  if (ff->type == FF_LENNARD || ff->type == FF_BUCKINGHAM)
    {
/* get appropriate type */
    ff_gromacs = gromacs_type(ff);

/* print LJ/Buck */
    fprintf(fp, "%s %s %d", ff->atom[0], ff->atom[1], ff_gromacs->type);
    for (i=0 ; i<ff_gromacs->data_expected ; i++)
      fprintf(fp, " %f", ff_gromacs->data[i]);
    fprintf(fp, "\n");

    }
  }

return(0);
}

/***********************************************************/
/* topology ie molecule, connectivity, force fields (.top) */
/***********************************************************/
gint write_gromacs_top(gchar *filename, struct model_pak *model)
{
FILE *fp;

fp = fopen(filename, "wt");

/* for each molecule write topology entry */
/* TODO - sort into unique molecules (atoms must also be ordered the same) */

write_gromacs_header(fp, model);

write_gromacs_atoms(fp, model);

write_gromacs_bonds(fp, model);

write_gromacs_angles(fp, model);

write_gromacs_dihedrals(fp, model);

write_gromacs_nonbond(fp, model);

write_gromacs_footer(fp, model);

fclose(fp);

return(0);
}

/****************/
/* file writing */
/****************/
gint write_gromacs(gchar *filename, struct model_pak *model)
{
gchar *temp;

temp = parse_extension_set(filename, "top");

write_gromacs_gro(filename, model);
write_gromacs_top(temp, model);

g_free(temp);

return(0);
}

/********************************/
/* read in a GROMACS FF library */
/********************************/
#define DEBUG_GROMACS_READ_FF 0
GSList *gromacs_read_ff(const gchar *filename)
{
gint state, type, num_tokens, ub, ua, ud, i;
gchar *line, **buff, **abuff, **dbuff;
GSList *list=NULL;
struct forcefield_pak *ff;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(NULL);

/* unknown bond, angle, dihedral types */
ub = ua = ud = 0;

/* stop auto skip of lines starting with # ... ie #define */
file_skip_comment(FALSE);

line = file_read_line(fp);
state = -1;
while (line)
  {
/* TODO - get rid of leading whitespace, so funny comment lines are properly ignored */

  ff = NULL;

/* directive processing */
  switch(line[0])
    {
    case '[':
      state = -1;
      if (g_strrstr(line, "bond"))
        state = BOND;
      if (g_strrstr(line, "angle"))
        state = ANGLE;
      if (g_strrstr(line, "dihedral"))
        state = DIHEDRAL;
      break;

/* process #defines */
    case '#':
      if (g_ascii_strncasecmp("#define", line, 7) == 0) 
        {
        buff = get_tokens(line, 3);
        abuff = g_strsplit(*(buff+1), "_", -1);

#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 6
        num_tokens = g_strv_length(abuff);
#else
        num_tokens = 0;
#endif

        if (abuff)
          {
          if (g_ascii_strncasecmp(*(abuff), "improper", 8) == 0)
            {
            if (num_tokens > 4)
              {
              ff = ff_type_new(FF_IMPROPER);
              g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(abuff+num_tokens-4));
              g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(abuff+num_tokens-3));
              g_snprintf(ff->atom[2], FF_MAX_SYMBOL, "%s", *(abuff+num_tokens-2));
              g_snprintf(ff->atom[3], FF_MAX_SYMBOL, "%s", *(abuff+num_tokens-1));
/* TODO - Y, Z -> X ... since X is wildcard type */
              for (i=4 ; i-- ; )
                if (ff->atom[i][0] == 'Y' || ff->atom[i][0] == 'Z')
                  ff->atom[i][0] = 'X';

              dbuff = get_tokens(line, 6);
              ff->bond_value = str_to_float(*(dbuff+2));
              ff->bond_units = FF_DEG;
              ff->data[0] = str_to_float(*(dbuff+3));
              ff->data[1] = str_to_float(*(dbuff+4));
              g_strfreev(dbuff);

              ff->atoms_current = ff->atoms_expected;
              ff->data_current = ff->data_expected = 2;

              list = g_slist_prepend(list, ff);
              ff = NULL;
              }
            }
          g_strfreev(abuff);
          }
        g_strfreev(buff);
        }
      state = -1;
      break;

/* comment */
    case ';':
      break;

    default:
/* normal processing */
      buff = tokenize(line, &num_tokens);
      switch (state)
        {
        case BOND:
          if (num_tokens > 4)
            {
            type = g_ascii_strtod(*(buff+2), NULL);
            switch (type)
              {
              case 1:
                ff = ff_type_new(FF_HARMONIC);
                g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(buff+0));
                g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(buff+1));
/* NB: GROMACS uses nm rather than angs */
                ff->bond_value = 10.0 * str_to_float(*(buff+3));
                ff->bond_units = FF_ANG;
/* nm -> ang correction */
                ff->data[0] = 0.01 * str_to_float(*(buff+4));
                ff->data_units = FF_KJ;

                ff->atoms_current = ff->atoms_expected;
                ff->data_current = ff->data_expected;
                break;

              case 3:
g_assert(num_tokens > 5);
                ff = ff_type_new(FF_MORSE);
                g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(buff+0));
                g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(buff+1));
/* nm -> ang correction */
                ff->bond_value = 10.0*str_to_float(*(buff+3));
                ff->bond_units = FF_ANG;
                ff->data[0] = str_to_float(*(buff+4));
                ff->data[1] = 0.1 * str_to_float(*(buff+5));
                ff->data_units = FF_KJ;

                ff->atoms_current = ff->atoms_expected;
                ff->data_current = ff->data_expected;
                break;

              default:
                ub++;
              }
            }
          break;

        case ANGLE:
          if (num_tokens > 5)
            {
            type = g_ascii_strtod(*(buff+3), NULL);
            switch (type)
              {
              case 1:
                ff = ff_type_new(FF_3B_HARMONIC);
                g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(buff+0));
                g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(buff+1));
                g_snprintf(ff->atom[2], FF_MAX_SYMBOL, "%s", *(buff+2));
                ff->bond_value = str_to_float(*(buff+4));
                ff->bond_units = FF_DEG;
                ff->data[0] = str_to_float(*(buff+5));
/* strictly - kJ/mol rad-2 */
                ff->data_units = FF_KJ;
                ff->atoms_current = ff->atoms_expected;
                ff->data_current = ff->data_expected;
                break;

              default:
                ua++;
              }
            }
          break;

        case DIHEDRAL:
          if (num_tokens > 4)
            {
            type = g_ascii_strtod(*(buff+4), NULL);
            switch (type)
              {
              case 1:
                if (num_tokens < 7)
                  {
                  ud++;
                  break;
                  }
                ff = ff_type_new(FF_DIHEDRAL);
                g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(buff+0));
                g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(buff+1));
                g_snprintf(ff->atom[2], FF_MAX_SYMBOL, "%s", *(buff+2));
                g_snprintf(ff->atom[3], FF_MAX_SYMBOL, "%s", *(buff+3));
                ff->bond_value = str_to_float(*(buff+5));
                ff->bond_units = FF_DEG;
                ff->data[0] = str_to_float(*(buff+6));
                ff->data_units = FF_KJ;
                ff->atoms_current = ff->atoms_expected;
                ff->data_current = ff->data_expected;
                break;

              case 3:
                if (num_tokens < 11)
                  {
                  ud++;
                  break;
                  }
                ff = ff_type_new(FF_DIHEDRAL_RB);
                g_snprintf(ff->atom[0], FF_MAX_SYMBOL, "%s", *(buff+0));
                g_snprintf(ff->atom[1], FF_MAX_SYMBOL, "%s", *(buff+1));
                g_snprintf(ff->atom[2], FF_MAX_SYMBOL, "%s", *(buff+2));
                g_snprintf(ff->atom[3], FF_MAX_SYMBOL, "%s", *(buff+3));
                ff->data[0] = str_to_float(*(buff+5));
                ff->data[1] = str_to_float(*(buff+6));
                ff->data[2] = str_to_float(*(buff+7));
                ff->data[3] = str_to_float(*(buff+8));
                ff->data[4] = str_to_float(*(buff+9));
                ff->data[5] = str_to_float(*(buff+10));
                ff->data_units = FF_KJ;

                ff->atoms_current = ff->atoms_expected;
                ff->data_current = ff->data_expected;
                break;

              default:
                ud++;
              }
            }
          break;
        }

/* add the object to the return list */
      if (ff)
        list = g_slist_prepend(list, ff);

      g_strfreev(buff);
      break;
    }

  g_free(line);
  line = file_read_line(fp);
  }

g_free(line);

#if DEBUG_GROMACS_READ_FF
printf("processed: %d items.\n", g_slist_length(list));
printf("ignored: [bonds = %d] [angles = %d] [dihedrals = %d]\n", ub, ua, ud);
/*
ff_dump_type(FF_IMPROPER, list);
*/
#endif

file_skip_comment(TRUE);

return(list);
}
