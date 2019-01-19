/*
 Copyright (C) 2005 by Menno Deij
 
 m.deij@science.ru.nl
 
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
#include <ctype.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "scan.h"
#include "matrix.h"
#include "interface.h"
#include "numeric.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

#define DEBUG_TOKENIZE_BONDLINE 0

gchar **tok_this_line(const gchar *src, gint *num)
{
  gint i, j, n, len;
  gchar *tmp, *ptr;
  gchar **dest;
  gchar *my_delim;
  GSList *list=NULL, *item=NULL;
  
  /* checks */
  if (!src)
    return(NULL);
  
  my_delim = g_strdup("[]");
  
  /* duplicate & replace all whitespace with a space */
  tmp = g_strdup(src);
  for (i=0 ; i<strlen(tmp) ; i++)
    if (isspace((int) *(tmp+i)) || *(tmp+i) == my_delim[0] || *(tmp+i) == my_delim[1])
      *(tmp+i) = ' ';
  
  /* strange errors can be avoided if a strstrip is done */
  g_strstrip(tmp);
  
#if DEBUG_TOKENIZE_BONDLINE
  printf("tokenizing [%s]\n", tmp);
#endif
  
  len = strlen(tmp);
  i=n=0;
  while(i<len)
  {
    /* find end of current token */
    j=i;
    while(!isspace((int) *(tmp+j)) && j<len)
      j++;
    
    /* assign token */
    ptr = g_strndup(tmp+i, j-i);
    
    list = g_slist_prepend(list, ptr);
    n++;
    
    /* find start of new token */
    i=j;
    
    while(isspace((int) *(tmp+i)) && i<len)
      i++;
  }
  list = g_slist_reverse(list);
  
  /* return a NULL if no tokens were found */
  if (!n)
  {
    *num = 0;
    g_free(tmp);
    free_slist(list);
    return(NULL);
  }
  
  /* num+1 -> last ptr is NULL, so g_strfreev works */
  dest = g_malloc((n+1)*sizeof(gchar *));
  
  i=0;
  /* fill in the non empty tokens */
  item = list;
  while (i<n)
  {
    if (item != NULL)
    {
      /* comment character - ignore all subsequent tokens */
      ptr = item->data;
      if (*ptr == '#')
        break;
      
      *(dest+i) = g_strdup(ptr);
#if DEBUG_TOKENIZE_BONDLINE
      printf(" (%s)", ptr);
#endif
      item = g_slist_next(item);
    }
    else
    {
      /* fake item */
      *(dest+i) = g_strdup(" ");;
#if DEBUG_TOKENIZE_BONDLINE
      printf(" (empty token)");
#endif
    }
    i++;
  }
  
  /* terminate */
  *(dest+i) = NULL;
  *num = i;
  
#if DEBUG_TOKENIZE_BONDLINE
  printf(" %p",*(dest+i));
  printf(": found %d tokens\n", *num);
#endif
  
  /* done */
  g_free(my_delim);
  g_free(tmp);
  free_slist(list);
  
  return(dest);
}

void frac_to_cartesian(gpointer core_ptr, gpointer data_ptr)
{
  struct core_pak *core = core_ptr;
  struct model_pak *data = data_ptr;
  
  vecmat(data->latmat, core->x);
}

void print_core_coordinates(gpointer core_ptr, gpointer dummy)
{
  struct core_pak *core = core_ptr;
  g_message("[ %3.2f %3.2f %3.2f ]", core->x[0], core->x[1], core->x[2]);
}

gchar **tokenize_bondline(FILE *fp, gint *num_tokens)
{
  gchar **buff, line[LINELEN];

  do
  {
    if (fgetline(fp, line))
    {
      *num_tokens = 0;
      return(NULL);
    }
    buff = tok_this_line(line, num_tokens);
  }
  while (!buff);
  
  return(buff);
}

/****************/
/* file writing */
/****************/
gint write_cgf(gchar *filename, struct model_pak *data)
{
  GSList *core_iter;/*FIX 7d2047*/
  GSList *primary_cores = NULL;
  FILE *fp;
  gint nr_primary_cores = 0, max_bonds = 0;/*FIX 7d2047*/
#ifdef UNUSED_BUT_SET
  GSList *bond_iter;/*FIX 7d2047*/
  gint bond_counter = 0;
  struct bond_pak *bond;
  gint offset_counter = 0;
#endif
   
  gint *nr_bonds = g_new(gint, g_slist_length(data->cores));
  
  struct core_pak *core;
  
  /* checks and open file */
  g_return_val_if_fail(data != NULL, 1);
  fp = fopen(filename, "wt");
  g_return_val_if_fail(fp != NULL, 2);
  
  fprintf(fp, "--------------------------------------------\n");
  fprintf(fp, "Title: CGF created by GDIS v%3.2f\n", VERSION);
  fprintf(fp, "a= %3.5f b= %3.5f c= %3.5f\n", data->pbc[0], data->pbc[1], data->pbc[2]);
  fprintf(fp, "alpha= %3.5f beta= %3.5f gamma= %3.5f\n", data->pbc[3] * R2D, data->pbc[4] * R2D, data->pbc[5] * R2D);
  fprintf(fp, "Spacegroup information:      SPGR = %-24s OPT = %i\n", data->sginfo.spacename, data->sginfo.cellchoice);
  
  gint total_bonds = 0;
  
  for (core_iter = data->cores; core_iter; core_iter = core_iter->next)
  {
    core = (struct core_pak *)core_iter->data;
    
    if ( (core->x[0]  < 1.0) && 
         (core->x[0] >= 0.0) &&
         (core->x[1]  < 1.0) && 
         (core->x[1] >= 0.0) &&
         (core->x[2]  < 1.0) && 
         (core->x[2] >= 0.0) )
    {
      nr_bonds[nr_primary_cores] = g_slist_length(core->bonds);
      total_bonds += nr_bonds[nr_primary_cores];
      max_bonds = nr_bonds[nr_primary_cores] > max_bonds ? nr_bonds[nr_primary_cores] : max_bonds;
      primary_cores = g_slist_prepend(primary_cores, core);
      ++nr_primary_cores;
    }
  }
  
#ifdef UNUSED_BUT_SET
  gdouble * bond_strengths = g_new0(gdouble, total_bonds);
  gint * bond_to_values =  g_new0(gint, total_bonds);
  gint * bond_offsets = g_new0(gint, total_bonds * 3);
#endif /*FIX 7d2047*/
  
  fprintf(fp, "Nr of centres of mass: %i\n", nr_primary_cores);
  fprintf(fp, "Maximum connectivity: %i\n", max_bonds);
  fprintf(fp, "These are fractional coordinates\n");
  fprintf(fp, "--------------------------------------------\n");

#ifdef UNUSED_BUT_SET
  bond_counter = 0;/*FIX 7d2047*/
  offset_counter = 0;

  
  for (core_iter = data->cores; core_iter; core_iter = core_iter->next)
  {
    core = core_iter->data;
    for (bond_iter = core->bonds; bond_iter; bond_iter = bond_iter->next)
    {
      bond = bond_iter->data;
      bond_to_values[bond_counter] = g_ascii_strtoull(g_strndup(bond->atom2->atom_label, 2), NULL, 10);
      bond_counter++;
    }
  }
//  gint i = 0;
//  for (; i < total_bonds; ++i)
//    g_message("%3.5f", bond_strengths[i]);
#endif /*FIX 7d2047*/  
  g_free(nr_bonds);
#ifdef UNUSED_BUT_SET
  g_free(bond_strengths);
  g_free(bond_to_values);
  g_free(bond_offsets);
#endif /*FIX 7d2047*/

  fclose(fp);
  return(0);
}

gint read_cgf(gchar *filename, struct model_pak *data)
{
  GSList *core_iter;
  gint n=0, num_tokens, num_gu, max_bonds, nr_bond_lines, current_bond_line;
  gint buff_pointer_offset, current_bond, w;
  gchar **buff, line[LINELEN];
  gchar *spacename = NULL;

  struct core_pak *primary_bondto_core, *sym_bondto_core, *temp_core;//, *found_core;
  struct core_pak *core;

  FILE *fp;
  
  gint *bondto; 
  gdouble *bond_offsets;
  
  /* checks */
  g_return_val_if_fail(data != NULL, 1);
  g_return_val_if_fail(filename != NULL, 1);
  
  fp = fopen(filename, "rt");
  if (!fp)
    return(1);
  
  fgetline(fp, line);
  g_strstrip(line);
  if (g_ascii_strcasecmp("--------------------------------------------", line) != 0)
  {
    /* not a valid crystal graph file*/
    gui_text_show(ERROR, "Incorrect format for Crystal Graph (.cgf) file\n");
    return(1);
  }
  
  /* crystal graphs are always 3d periodic */
  data->periodic = 0;
  /* particle locations should always be fractional  *
  data->fractional = FALSE;
  */
  fgetline(fp, line);

  /* no problems with char * and gchar * ?? */
  /* the next line give a bus error  */
  //sscanf(line, "Title: %20[a-zA-Z0-9/ ]", data->title); 

  /* next line contains created by information, not used here */
  fgetline(fp, line);

  /* next line contains axis length information as a, b, c */
  buff = get_tokenized_line(fp, &num_tokens);
  data->pbc[0] = str_to_float(*(buff+1));
  data->pbc[1] = str_to_float(*(buff+3));
  data->pbc[2] = str_to_float(*(buff+5));
  g_strfreev(buff);

  /* next line contains axis angle information as alpha, beta, gamma*/
  buff = get_tokenized_line(fp, &num_tokens);
  data->pbc[3] = str_to_float(*(buff+1));
  data->pbc[4] = str_to_float(*(buff+3));
  data->pbc[5] = str_to_float(*(buff+5));
  g_strfreev(buff);

  /* convert angles to radians  */
  data->pbc[3] *= D2R; data->pbc[4] *= D2R; data->pbc[5] *= D2R;

  /* create cartesian axes */
  matrix_lattice_init(data);

  /* spacegroup info */
  fgetline(fp, line); 
  sscanf(line, "Spacegroup information:      SPGR = %11[a-zA-Z0-9-/ ] OPT =%2d",
         spacename, &data->sginfo.cellchoice);

  g_strstrip(spacename);
  data->sginfo.spacename = g_strdup(spacename);
  
  /* next line contains number of GUs in the unit cell  */
  fgetline(fp, line);
  sscanf(line, "Nr of centres of mass: %2d", &num_gu);

#ifdef DEBUG_CGF_READ
  g_message("%d number of mass", num_gu);
#endif
  
  fgetline(fp, line);
  sscanf(line, "Maximum connectivity: %2d", &max_bonds);

#ifdef DEBUG_CGF_READ
  g_message("%d bonds", max_bonds);
#endif

  /* for each 5 bonds there is a line of bond specifications  */
  nr_bond_lines = max_bonds / 5 + ( max_bonds % 5 == 0 ? 0 : 1);
  
#ifdef DEBUG_CGF_READ
  g_message("%d lines of bond information", nr_bond_lines);
#endif
  
  bondto = g_new0(gint, nr_bond_lines * num_gu * 5);/*FIX e10cdb*/
  bond_offsets = g_new(gdouble, nr_bond_lines * num_gu * 15);
  //gdouble *bond_strengths;
  //bond_strengths = g_new(gdouble, nr_bond_lines * num_gu * 5);
  
  gint bondto_counter = 0, bond_offsets_counter = 0;
  
  /* create GUs a priori  */
  for (n = 0; n < num_gu; ++n)
  {
    core = core_new("C", NULL, data); 
    data->cores = g_slist_append(data->cores, core);
    if (core->atom_label)
      g_free(core->atom_label);
    core->atom_label = g_strdup_printf("%i [0 0 0]", n);
  }
  
  /* skip the next 2 lines  */
  for (n = 0; n < 2; ++n)
    fgetline(fp,line);
 
  /* TODO: for now, only GUs are read in, bonds will be done later  */
  for (n = 0; n < num_gu; ++n)
  {
    /* get the nth core from the list  */
    core = g_slist_nth_data(data->cores, n);

    for (current_bond_line = 0; current_bond_line < nr_bond_lines; ++current_bond_line)
    {      
      /* first line contains bond-to info  */
      buff = get_tokenized_line(fp, &num_tokens);
      
      for (current_bond = 0; current_bond < 5; ++current_bond)
      {
        gdouble to = g_ascii_strtod(*(buff + current_bond), NULL);
        bondto[bondto_counter++] = nearest_int(to);
      }
      
      g_strfreev(buff);
      
      /* second line contains GU fractional coordinates and bond offsets  */
      buff = tokenize_bondline(fp, &num_tokens);
      
      buff_pointer_offset = 0;
      if (current_bond_line == 0)
      {          
        core->x[0] = str_to_float(*(buff+2));
        core->x[1] = str_to_float(*(buff+3));
        core->x[2] = str_to_float(*(buff+4));
        
        /* in this case the buff tokenizer contains 3 more tokens */
        buff_pointer_offset = 5;
      }
      

      for (current_bond = 0; current_bond < 15; current_bond++)
      {
        bond_offsets[bond_offsets_counter++] = g_ascii_strtod(*(buff + buff_pointer_offset + current_bond), NULL);
      }
      
      g_strfreev(buff);
      
      /* the next line contains the labels of the bonds  */
      fgetline(fp, line);
      
      /* the next line contains the bond strengths *
      buff = tokenize_bondline(fp, &num_tokens);
      for (current_bond = 0; current_bond < 5; current_bond++)
      {
        bond_strengths[bond_strength_counter++] = g_ascii_strtod(*(buff+ current_bond), NULL);
      }
      g_strfreev(buff); 
      */
      /* the next line is blank  */
      fgetline(fp, line);
      
    }
    /* between particle specifications one extra blank line is present */
    fgetline(fp, line);
  }

  gint total_bonds = 0;
  /* now we create new cores for each bond */
  w = 0;
  for (n = 0; n < num_gu; ++n)
  {
    /* get the core that the bonds are coming from */
    core = g_slist_nth_data(data->cores, n);    
    
    for (current_bond = 0; current_bond < nr_bond_lines * 5; ++current_bond)
    {
      /* compute the correct index for the bondto array */
      gint bondto_index = current_bond + n * nr_bond_lines * 5;
      
      if (bondto[bondto_index] != 0)
      {
        total_bonds++;
        /* get the core that has the identity of the core that the bond goes to */
        primary_bondto_core = g_slist_nth_data(data->cores, bondto[bondto_index]-1);

        /* if that core is inside the unit cell, do not create a new one */
        if (bond_offsets[w] == 0 && bond_offsets[w+1] == 0 && bond_offsets[w+2] ==0)
        {
          sym_bondto_core = primary_bondto_core;
          w+=3;
        }
        else 
        {
          gchar * primary_label = primary_bondto_core->atom_label;
          gchar * sym_label = g_strdup_printf("%s [%i %i %i]", 
                                               g_strndup(primary_label, 1), 
                                               (int)bond_offsets[w], 
                                               (int)bond_offsets[w+1], 
                                               (int)bond_offsets[w+2]);
          
          sym_bondto_core = NULL;
          for (core_iter = data->cores; core_iter; core_iter = core_iter->next)
          {
            temp_core = core_iter->data;
            if (g_ascii_strcasecmp(temp_core->atom_label, sym_label) == 0)
            { 
              sym_bondto_core = temp_core;
            }
          }
          
          if (sym_bondto_core == NULL)
          {
            sym_bondto_core = core_new("C", NULL, data);
            sym_bondto_core->atom_label = sym_label;
            data->cores = g_slist_append(data->cores, sym_bondto_core);
          }

          /* set the new core's coordinates to the core it represents
            and add the offset of the bond */
          sym_bondto_core->x[0] = primary_bondto_core->x[0] + bond_offsets[w++];          
          sym_bondto_core->x[1] = primary_bondto_core->x[1] + bond_offsets[w++];          
          sym_bondto_core->x[2] = primary_bondto_core->x[2] + bond_offsets[w++];
        }
        /* now, connect the bond between the core that the bond comes
          from and the (possibly) newly created core that the bond goes to */
        connect_user_bond(core, sym_bondto_core, BOND_SINGLE, data);
      }
      /* the  bond goes to 0, which means that there is no bond */
      else 
      {  
        w += 3;
      }      
    }
  }

  /* now that we have all cores, their coordinates need to be transformed to cartesian */
  g_slist_foreach(data->cores, frac_to_cartesian, data);
  
  g_free(bondto);
  g_free(bond_offsets);
  
  /* got everything */
  data->num_asym = data->num_atoms = num_gu;
  
  /* model setup */
  strcpy(data->filename, filename); 
  g_free(data->basename);
  data->basename = parse_strip(filename);

  model_prep(data);

  fclose(fp);
  return (0);
}
