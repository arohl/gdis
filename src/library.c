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
#include "library.h"
#include "matrix.h"
#include "model.h"
#include "scan.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/*************************/
/* free a library folder */
/************************/
void library_folder_free(gpointer data)
{
GSList *list;
struct folder_pak *folder = data;
struct entry_pak *entry;

list = folder->list;
while (list)
  {
  entry = list->data;
  list = g_slist_next(list);
  g_free(entry->name);
  g_string_free(entry->info, TRUE);
  g_free(entry->offset);
  g_free(entry);
  }
g_slist_free(folder->list);
g_free(folder);
}

/******************************/ 
/* create a new library entry */
/******************************/ 
gpointer library_entry_add(gchar *folder_name, gchar *entry_name, gpointer offset)
{
struct folder_pak *folder;
struct entry_pak *entry;

/* TODO - check for entry repetitions within a given folder */
entry = g_malloc(sizeof(struct entry_pak));
entry->offset = offset;
entry->name = g_strdup(entry_name);
entry->info = g_string_new(NULL);

folder = g_hash_table_lookup(sysenv.library, folder_name);
if (!folder)
  {
  folder = g_malloc(sizeof(struct folder_pak));
  folder->list = NULL;
  g_hash_table_insert(sysenv.library, g_strdup(folder_name), folder);
  }
folder->list = g_slist_prepend(folder->list, entry);

return(entry);
}

/********************************/
/* setup the library hash table */
/********************************/
void library_init(void)
{
gint num_tokens;
gchar **buff, *folder=NULL, *filename;
fpos_t *offset;
struct entry_pak *entry = NULL;
FILE *fp;

/* setup hash table */
sysenv.library = g_hash_table_new_full(g_str_hash, g_str_equal,
                                       g_free, library_folder_free);

/* scan the library file (gdis.library = cif format) for folder/entry info */
filename = g_build_filename(sysenv.gdis_path, LIBRARY_FILE, NULL);

printf("scanning: %s\n", filename);

fp = fopen(filename, "rt");
g_free(filename);
if (!fp)
  return;
for (;;)
  {
  buff = get_tokenized_line(fp, &num_tokens);
  if (!num_tokens)
    break;
  if (num_tokens < 2)
    {
    g_strfreev(buff);
    continue;
    }

/* keyword search */
  if (g_ascii_strncasecmp(*buff, "_gdis_folder_name", 17) == 0)
    {
    folder = g_strdup(*(buff+1));
    entry = NULL;
    }

  if (g_ascii_strncasecmp(*buff, "_gdis_entry_name", 16) == 0)
    {
    offset = g_malloc(sizeof(fpos_t));
    fgetpos(fp, offset);

    if (folder)
      {
      entry = library_entry_add(folder, *(buff+1), offset);
      g_free(folder);folder=NULL;/*FIX 3fe0a9*/
      }
    else
      entry = library_entry_add("default", *(buff+1), offset);
    }

  if (g_ascii_strncasecmp(*buff, "_database_code_ICSD", 19) == 0)
    {
    if (entry)
      g_string_append_printf(entry->info, "ICSD code %s\n", *(buff+1));
    }

  if (g_ascii_strncasecmp(*buff, "_journal_year", 13) == 0)
    {
    if (entry)
      g_string_append_printf(entry->info, "Published in %s\n", *(buff+1));
    }

  if (g_ascii_strncasecmp(*buff, "_refine_ls_R_factor", 19) == 0)
    {
    if (entry)
      g_string_append_printf(entry->info, "Least squares R factor %s\n", *(buff+1));
    }

  g_strfreev(buff);
  }
fclose(fp);

/* TODO - foreach folder, sort the entry list alphabetically? */

}

/***************************/
/* retrive a library entry */
/***************************/
#define DEBUG_LOAD_CIF 0
gint library_entry_get(gpointer offset, struct model_pak *data)
{
gint i, j, n, first, new, order, pos, keyword, loop_count=0;
gint min_tokens, num_tokens, len, flag;
gint s_pos, l_pos, x_pos, y_pos, z_pos, o_pos;
gchar **buff, *filename, *tmp, *name=NULL, line[LINELEN];
gdouble sign, mat[9], off[3];
GSList *item, *list=NULL;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);

filename = g_build_filename(sysenv.gdis_path, LIBRARY_FILE, NULL);
fp = fopen(filename, "rt");
g_free(filename);
if (!fp)
  return(2);
if (fsetpos(fp, offset))
  {
  printf("Error positioning file pointer.\n");
  return(3);
  }

/* atom counter */
n=-1;
new = first = 0;
data->id = -1;
for(;;)
  {
/* end if EOF */
  if (fgetline(fp, line))
    break;

/* search for data */
  list = get_keywords(line);
cif_parse:;
  if (list != NULL)
    {
    keyword = GPOINTER_TO_INT(list->data);
    switch(keyword)
      {
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
/* LIBRARY - terminate after 1 model loaded */
goto cif_done;

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
        data->pbc[3] = D2R*str_to_float(*(buff+1));
        g_strfreev(buff);
        break;
      case CIF_CELL_BETA:
        buff = get_tokens(line, 3);
        data->pbc[4] = D2R*str_to_float(*(buff+1));
        g_strfreev(buff);
        break;
      case CIF_CELL_GAMMA:
        buff = get_tokens(line, 3);
        data->pbc[5] = D2R*str_to_float(*(buff+1));
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
          if (fgetline(fp, line))
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
                  strncpy(line, *(buff+i+n)+j-1, 3);
                  line[3] = '\0';
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
/* while keyword lists are found */
        list = get_keywords_anywhere(line);
        while (list != NULL)
          {
/* enumerate keyword list */
          for (item=list ; item ; item=g_slist_next(item))
            {
            keyword = GPOINTER_TO_INT(item->data);
            switch(keyword)
              {
              case CIF_TYPE_SYMBOL:
                s_pos = loop_count;
                break;
              case CIF_LABEL:
                l_pos = loop_count;
                break;
              case CIF_CART_X:
                data->fractional = FALSE;
                x_pos = loop_count;
                break;
              case CIF_CART_Y:
                data->fractional = FALSE;
                y_pos = loop_count;
                break;
              case CIF_CART_Z:
                data->fractional = FALSE;
                z_pos = loop_count;
                break;
              case CIF_FRAC_X:
                data->fractional = TRUE;
                x_pos = loop_count;
                break;
              case CIF_FRAC_Y:
                data->fractional = TRUE;
                y_pos = loop_count;
                break;
              case CIF_FRAC_Z:
                data->fractional = TRUE;
                z_pos = loop_count;
                break;
              case CIF_SOF:
                o_pos = loop_count;
                break;
              }
            }

/* get next line and keyword list */
          loop_count++;
          if (fgetline(fp, line))
            goto cif_done;
          g_slist_free(list); 
          list = get_keywords_anywhere(line);
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

/* NEW - not enough data for valid atoms => ignore this loop_ (but don't quit) */
/*
          goto cif_done;
*/
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
          if (fgetline(fp, line))
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

