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
#include <string.h>
#include <time.h>

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "model.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***************/
/* CIF writing */
/***************/
gint write_cif(gchar *filename, struct model_pak *data)
{
gint flag=0;
gdouble tmat[9], x[3], depth=1.0;
GSList *list;
FILE *fp;
time_t t1;
struct core_pak *core;

/* init & checks */
g_return_val_if_fail(data != NULL, 1);
fp = fopen(filename, "wt");
g_return_val_if_fail(fp != NULL, 2);
  
/* is it a surface with negative z? */
if (data->periodic == 2)
  {
  if (g_slist_length(data->cores))
    {
    core = (struct core_pak *) g_slist_nth_data(data->cores, 0);
    if (core->x[2] < 0.0)
      flag++;
    }
  }

/* rot to make z positive */
matrix_rotation(&tmat[0], PI, PITCH);

t1 = time(NULL);
fprintf(fp, "data_block_1\n");
fprintf(fp, "_audit_creation_date         '%s'\n", g_strstrip(ctime(&t1)));
fprintf(fp, "_audit_creation_method       'generated by GDIS v%4.2f'\n", VERSION);
fprintf(fp, "\n\n");
if (data->periodic)
  {
  fprintf(fp, "_cell_length_a            %8.4f\n",data->pbc[0]);
  fprintf(fp, "_cell_length_b            %8.4f\n",data->pbc[1]);

if (data->periodic == 2)
  {
/* get depth info */
  depth = (data->surface.region[0]+data->surface.region[1])*data->surface.depth;
/* no depth info - make it large enough to fit everything */
  if (depth < POSITION_TOLERANCE)
    depth = 2.0*data->rmax;
  fprintf(fp, "_cell_length_c            %8.4f\n", depth);
  }
else fprintf(fp, "_cell_length_c            %8.4f\n",data->pbc[2]);

  fprintf(fp, "_cell_angle_alpha       %8.2f\n",180.0*data->pbc[3]/PI);
  fprintf(fp, "_cell_angle_beta        %8.2f\n",180.0*data->pbc[4]/PI);
  fprintf(fp, "_cell_angle_gamma       %8.2f\n",180.0*data->pbc[5]/PI);
  fprintf(fp, "\n\n");
  fprintf(fp, "_symmetry_space_group_name_H-M   '%s'\n",data->sginfo.spacename);
  fprintf(fp, "_symmetry_Int_Tables_number       %d\n",data->sginfo.spacenum);
  fprintf(fp, "\n\n");
  }

/* coords - format section */
fprintf(fp, "loop_\n");
fprintf(fp, "_atom_site_label\n");
if (data->periodic)
  {
  fprintf(fp, "_atom_site_fract_x\n");
  fprintf(fp, "_atom_site_fract_y\n");
  fprintf(fp, "_atom_site_fract_z\n");
  }
else
  {
  fprintf(fp, "_atom_site_cartn_x\n");
  fprintf(fp, "_atom_site_cartn_y\n");
  fprintf(fp, "_atom_site_cartn_z\n");
  }

/* coords - data section */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* only the asymmetric cell if periodic */
  if (data->periodic && !core->primary)
    continue;

/* NB: want fractional if 3D periodic, otherwise cartesian */
  ARR3SET(x, core->x);
/* transformation needed? */
    if (flag)
      vecmat(tmat, x);
/* make the z coordinate "fractional" for a surface */
  if (data->periodic == 2)
    x[2] /= depth;

  fprintf(fp, "%2s %9.4f %9.4f %9.4f\n", 
          elements[core->atom_code].symbol, x[0], x[1], x[2]);
  }
fprintf(fp, "\n\n");

fclose(fp);
return(0);
}

/***************/
/* CIF loading */
/***************/
#define DEBUG_LOAD_CIF 0
gint read_cif(gchar *filename, struct model_pak *data)
{
gint i, j, n, first, new, order, pos, keyword, loop_count=0;
gint min_tokens, num_tokens, len, flag;
gint s_pos, l_pos, x_pos, y_pos, z_pos, o_pos;
gchar **buff, *tmp, *name=NULL, *line;
gdouble sign, mat[9], off[3];
GSList *list=NULL;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

fp = fopen(filename, "rt");
if (!fp)
  return(3);

/* atom counter */
n=-1;
new = first = 0;
data->id = -1;
for(;;)
  {
/* end if EOF */
  line = file_read_line(fp);
  if (!line)
    break;

/*
  if (fgetline(fp, line))
    break;
*/

/* search for data */
  list = get_keywords(line);
cif_parse:;
  if (list != NULL)
    {
    keyword = GPOINTER_TO_INT(list->data);
    switch(keyword)
      {
/* model labels */
/* FIXME - needs to search for the 1st occurrance of ' or " & then get the string */
/*
      case CIF_CHEMICAL_NAME:
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++)
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
        g_free(data->basename);
        data->basename = g_strdup(g_strstrip(tmp));
        g_free(tmp);
        break;
      case CIF_MINERAL_NAME:
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++)
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
        g_free(data->basename);
        data->basename = g_strdup(g_strstrip(tmp));
        g_free(tmp);
        break;
*/

/* stopgap model name */
      case CIF_DATA_START:
        name = g_strdup(g_strstrip(&line[5]));
        break;

/* candidate new model trigger */
/* CIF is a pain - there seems to be no clear cut new structure */ 
/* trigger (in a multi structure file) that everyone uses */
      case CIF_AUDIT:
/* skip this the 1st time (we've already alloc'd a data pointer for safety) */
        new++;
        if ((new-first) > 1)
          {
/* NEW - some dodgy models have no atom data - so allow new model */
/* creation if we have (the bare minimum) some new cell parameters */
          if (!data->periodic)
            {
/* no cell parameters found yet - don't create a new model */
            new--;
            break;
            }
#  if DEBUG_LOAD_CIF
printf("Found %d atoms. [reset]\n", n);
#endif
/* alloc new pointer */
          data = model_new();
          if (data == NULL)
            goto cif_done;
          }

/* if we found a name - assign it now */
        if (name)
          {
          g_free(data->basename);
          data->basename = name;
          name = NULL;
          }
        
#if DEBUG_LOAD_CIF
printf("Start of new model: %s.\n", data->basename);
#endif
        break;

/* model specific data - the *data ptr MUST be allocated */
      case CIF_CELL_A:
        buff = get_tokens(line, 3);
        data->pbc[0] = str_to_float(*(buff+1));
        data->periodic++;
        g_strfreev(buff);
        break;
      case CIF_CELL_B:
        buff = get_tokens(line, 3);
        data->pbc[1] = str_to_float(*(buff+1));
        data->periodic++;
        g_strfreev(buff);
        break;
      case CIF_CELL_C:
        buff = get_tokens(line, 3);
        data->pbc[2] = str_to_float(*(buff+1));
        data->periodic++;
        g_strfreev(buff);
        break;
      case CIF_CELL_ALPHA:
        buff = get_tokens(line, 3);
        data->pbc[3] = PI*str_to_float(*(buff+1))/180.0;
        g_strfreev(buff);
        break;
      case CIF_CELL_BETA:
        buff = get_tokens(line, 3);
        data->pbc[4] = PI*str_to_float(*(buff+1))/180.0;
        g_strfreev(buff);
        break;
      case CIF_CELL_GAMMA:
        buff = get_tokens(line, 3);
        data->pbc[5] = PI*str_to_float(*(buff+1))/180.0;
        g_strfreev(buff);
        break;
      case CIF_SPACE_NAME:
/* remove the enclosing ' characters */
        tmp = g_strdup(get_token_pos(line,1));
        for (i=0 ; i<strlen(tmp) ; i++)
          if (*(tmp+i) == '\'')
            *(tmp+i) = ' ';
/* store the name, stripping spaces */
        data->sginfo.spacename = g_strdup(g_strstrip(tmp));
/* indicate that name should used in lookup */
        data->sginfo.spacenum = -1;
#if DEBUG_LOAD_CIF
printf("spacegroup: [%s]\n",data->sginfo.spacename);
#endif
        g_free(tmp);
        break;
      case CIF_SPACE_NUM:
        buff = get_tokens(line, 3);
        g_strfreev(buff);
        break;

      case CIF_EQUIV_SITE:
        loop_count++;
        break;

/* NEW - reinserted the lost symmetry matrix code */
      case CIF_EQUIV_POS:
/* allocate for order number of pointers (to matrices) */
        data->sginfo.matrix = (gdouble **) g_malloc(sizeof(gdouble *));
        data->sginfo.offset = (gdouble **) g_malloc(sizeof(gdouble *));
        data->sginfo.order=order=0;
        for(;;)
          {
/* terminate on EOF, */
          g_free(line);
          line = file_read_line(fp);
          if (!line)
            break;
/* blank line, */
          g_strstrip(line);
          if (!strlen(line))
            break;
/* or a new command */
          list = get_keywords(line);
          if (list)
            goto cif_parse;

/* TODO - make this parsing a subroutine */
          for(i=0 ; i<strlen(line) ; i++)
            if (*(line+i) == '\'' || *(line+i) == ',')
              *(line+i) = ' ';
          g_strstrip(line);
          buff = tokenize(line, &num_tokens);

          n = loop_count;
          while (n < num_tokens-2)
            {
/* FIXME - yet another mess that a linked list would greatly simplify */
/* number of ops */
          data->sginfo.matrix = (gdouble **) g_renew
                                (gdouble *, data->sginfo.matrix , (order+1));
          data->sginfo.offset = (gdouble **) g_renew
                                (gdouble *, data->sginfo.offset , (order+1));

/* actual op */
          *(data->sginfo.matrix+order) = (gdouble *) g_malloc(9*sizeof(gdouble));
          *(data->sginfo.offset+order) = (gdouble *) g_malloc(3*sizeof(gdouble));

#if DEBUG_LOAD_CIF
printf("[%s] [%s] [%s]\n", *(buff+n+0), *(buff+n+1), *(buff+n+2));
#endif

          VEC3SET(&mat[0], 0.0, 0.0, 0.0);
          VEC3SET(&mat[3], 0.0, 0.0, 0.0);
          VEC3SET(&mat[6], 0.0, 0.0, 0.0);
          VEC3SET(&off[0], 0.0, 0.0, 0.0);

          for (i=0 ; i<3 ; i++)
            {
            pos=0;
            sign=1.0;
            for (j=0 ; j<strlen(*(buff+i+n)) ; j++)
              {
              switch (*(*(buff+i+n)+j))
                {
                case '-':
                  sign = -1.0;
                  break;
                case '+':
                  sign = +1.0;
                  break;
                case 'x':
                  mat[i*3] = sign*1.0;
                  pos++;
                  break;
                case 'y':
                  mat[i*3 + 1] = sign*1.0;
                  pos++;
                  break;
                case 'z':
                  mat[i*3 + 2] = sign*1.0;
                  pos++;
                  break;
/* FIXME - a bit crude */
                case '/':

                  g_free(line);
                  line = g_strndup(*(buff+i+n)+j-1, 3);

                  off[i] = sign * str_to_float(line);
                  break;
                }
              }
            }
          ARR3SET((*(data->sginfo.matrix+order)+0), &mat[0]);
          ARR3SET((*(data->sginfo.matrix+order)+3), &mat[3]);
          ARR3SET((*(data->sginfo.matrix+order)+6), &mat[6]);
          ARR3SET((*(data->sginfo.offset+order)+0), &off[0]);

#if DEBUG_LOAD_CIF
P3MAT("output: ", *(data->sginfo.matrix+order));
P3VEC("output: ", *(data->sginfo.offset+order));
printf("\n\n");
#endif
          order++;
          data->sginfo.order++;
          n += 3;
          }
        g_strfreev(buff);
          }

#if DEBUG_LOAD_CIF
printf("Found %d symmetry matrices.\n", order);
#endif
        data->sginfo.order = order;
        break;

      case CIF_LOOP_START:
        loop_count = 0;
        break;

/* parsing for column info ie x,y,z positions frac/cart etc. */
      case CIF_ATOM_SITE:
        s_pos = l_pos = x_pos = y_pos = z_pos = o_pos = -1;

        while (g_strrstr(line, "atom_site") != NULL)
          {

          if (g_strrstr(line, "type_symbol"))
            s_pos = loop_count;
          if (g_strrstr(line, "label"))
            l_pos = loop_count;
          if (g_strrstr(line, "cart_x"))
            {
            x_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "cart_y"))
            {
            y_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "cart_z"))
            {
            z_pos = loop_count;
            data->fractional = FALSE;
            }
          if (g_strrstr(line, "fract_x"))
            {
            x_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "fract_y"))
            {
            y_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "fract_z"))
            {
            z_pos = loop_count;
            data->fractional = TRUE;
            }
          if (g_strrstr(line, "occupancy"))
            o_pos = loop_count;

/* get next line and keyword list */
          loop_count++;

          g_free(line);
          line = file_read_line(fp);

          if (!line)
            goto cif_done;
          }

/* either symbol or label can be present & used for atom id purposes */
        if (s_pos < 0)
          s_pos = l_pos;
        if (l_pos < 0)
          l_pos = s_pos;
/* check for minimum data */
        if (s_pos < 0 || x_pos < 0 || y_pos < 0 || z_pos < 0)
          {
#if DEBUG_LOAD_CIF
printf("read_cif() warning: incomplete cif file? [%d:%d:%d:%d]\n",
                                      s_pos, x_pos, y_pos, z_pos);
#endif
          break;
          }

/* the expected number of tokens */
        min_tokens = loop_count;

#if DEBUG_LOAD_CIF
printf(" min tokens: %d\n", min_tokens);
printf("data format: [%d] (%d) - [%d] [%d] [%d]  (%d)",
                s_pos, l_pos, x_pos, y_pos, z_pos, o_pos);
if (data->fractional)
  printf("(frac)\n");
else
  printf("(cart)\n");
#endif

/* while no new keywords found */
        n = 0;

        list = get_keywords(line);

        while (list == NULL)
          {
          buff = tokenize(line, &num_tokens);

/* NB: cif is stupid - it allows data to continue on the next line, */
/* until it gets its number of tokens */
/* hopefully, this'll allow us to get something, even on short lines */
          if (num_tokens >= z_pos)
            {

/* NEW - ignore labelled (*) symmetry equiv. atoms */
flag = TRUE;
if (l_pos >= 0 && l_pos < num_tokens)
  {
  len = strlen(*(buff+l_pos));
  if ((*(buff+l_pos))[len-1] == '*')
    flag = FALSE;
  }

            if (elem_symbol_test(*(buff+s_pos)) && flag)
              {
              if (l_pos >= 0 && l_pos < num_tokens)
                core = core_new(*(buff+s_pos), *(buff+l_pos), data);
              else
                core = core_new(*(buff+s_pos), NULL, data);

              data->cores = g_slist_prepend(data->cores, core);

              core->x[0] = str_to_float(*(buff+x_pos));
              core->x[1] = str_to_float(*(buff+y_pos));
              core->x[2] = str_to_float(*(buff+z_pos));

/* only get occupancy if we're sure we have enough tokens */
              if (o_pos > -1 && num_tokens >= min_tokens)
                {
                core->sof = str_to_float(*(buff+o_pos));
                core->has_sof = TRUE;
                }
              n++;
              }
            }
#if DEBUG_LOAD_CIF
          else
            printf("not enough tokens found.\n");
#endif
/* get next line */
          g_strfreev(buff);

          g_free(line);
          line = file_read_line(fp);
          if (!line)
            goto cif_done;

          list = get_keywords(line);
/* CURRENT - we really want this keyword list parsed again, just in case */
/* a new model trigger (ie audit_creation_date) was found */
          if (list)
            goto cif_parse;
          }
        break;
      } 
    }
  }
cif_done:;

g_free(line);

/* yet another hack to support the dodgy CIF standard */
j = new - first;
if (!j && n)
  {
/* no new model was triggered - but we found atoms, so we'll */
/* assume only one model was in the file & hope for the best */
  new++;
  }

#if DEBUG_LOAD_CIF
printf("Found %d atoms.\n", n);
printf("Found %d symmetry matrices.\n", data->sginfo.order);
printf("Found %d model(s)\n", new-first);
#endif

/* setup for display */
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  data = list->data;
  if (data->id == -1)
    {
    data->id = CIF;
    data->cores = g_slist_reverse(data->cores);
    model_prep(data);
    }
  }

#if DEBUG_LOAD_CIF
printf("Setup %d model(s).\n",new-first);
#endif

/* clean up & exit */
if (list)
  g_slist_free(list);

return(0);
}

