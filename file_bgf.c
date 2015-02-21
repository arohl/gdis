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

/*
 The BGF file type is the main file type for using GULP with the Dreiding force field.
 As OpenBabel can convert XYZ files to BGF files, this format is used to get both
 Dreiding atom types and Gasteiger charges.
 See also type.c for the typing routine 
 */

/****************/
/* file writing */
/****************/
gint write_bgf(gchar *filename, struct model_pak *data)
{
  GSList *core_iter, *bond_iter;
  FILE *fp;
  gint n = 0, nr_bonds = 0, bond_order;
  
  struct core_pak *core;
  struct bond_pak *bond;
  
  /* checks and open file */
  g_return_val_if_fail(data != NULL, 1);
  fp = fopen(filename, "wt");
  g_return_val_if_fail(fp != NULL, 2);
  
  fprintf(fp, "BIOGRF 200\n");
  fprintf(fp, "DESCRP %s\n", data->title);
  fprintf(fp, "REMARK BGF file created by GDIS v%3.2f\n", VERSION);
  fprintf(fp, "FORCEFIELD DREIDING   \n");
  fprintf(fp, "FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)\n");
  
  for (core_iter = data->cores; core_iter; core_iter = core_iter->next)
  {
    core = core_iter->data;
    /* replace any H_n GULP types back to H_ DREIDING type *
       NB. when the atom_type is not specified, this gives errors and all types are H_!*/
    if (core->atom_type == NULL)
      core->atom_type = g_strdup("UNK");
    if (g_ascii_strcasecmp(core->atom_type, "H_n") == 0)
    {
      g_free(core->atom_type);
      core->atom_type = g_strdup("H_"); 
    }
    fprintf(fp, "HETATM%6i %-6sRES A   444  % 3.4f  % 3.4f  % 3.4f %-4s   3 0  % 3.5f\n", 
            ++n, core->atom_label, core->x[0], core->x[1], core->x[2], core->atom_type, core->charge);
  }
  /*
  FORMAT CONECT (a6,12i6)
    
  CONECT     1     5    25     2
  ORDER      1     2     1     1
   */
  fprintf(fp, "FORMAT CONECT (a6,12i6)\n\n");
  for (n = 1, core_iter = data->cores; core_iter; core_iter = core_iter->next, ++n)
  {
    core = core_iter->data;
    /* n counts the number of bonds */
    nr_bonds = g_slist_length(core->bonds);
    if (nr_bonds > 0)
    {
      fprintf(fp, "CONECT%6i", n);
      for (bond_iter = core->bonds; bond_iter; bond_iter = bond_iter->next)
      {
        bond = bond_iter->data;
        
        fprintf(fp, "%6i", g_slist_index(data->cores, 
                                         (bond->atom2 == core) ? bond->atom1 : bond->atom2) + 1);
      }
      fprintf(fp, "\nORDER %6i", n);
      for (bond_iter = core->bonds; bond_iter; bond_iter = bond_iter->next)
      {
        bond = bond_iter->data;
        bond_order = 1;
        bond_order = (bond->type == BOND_DOUBLE) ? 2 : ( (bond->type == BOND_TRIPLE) ? 3 : 1);
        fprintf(fp, "%6i", bond_order);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "END\n");
  fclose(fp);
  return(0);
}

void species_line_from_hash_table(gpointer key, gpointer value, gpointer user_data)
{
  gchar *atom_type = (gchar *)key;
  gchar *atom_label = (gchar *)value;
  gchar **species_string = (gchar **)user_data;
  gchar *species_line = g_strdup_printf("%s 0.00000 %s", atom_label, atom_type);

  *species_string = g_strjoin("\n", *species_string, species_line, NULL);
  
  g_free(species_line);
}

/****************/
/* file reading */
/****************/

gint read_bgf(gchar *filename, struct model_pak *data)
{
  gint n=0, num_tokens;
  gchar **buff, line[LINELEN];
  struct core_pak *core;
  
  GHashTable * atom_types_and_labels = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
  
  FILE *fp;
  
  /* checks */
  g_return_val_if_fail(data != NULL, 1);
  g_return_val_if_fail(filename != NULL, 1);
  
  fp = fopen(filename, "rt");
  if (!fp)
    return(1);

/* old checks ... really necessary??? */
/*
  if (g_ascii_strcasecmp("BIOGRF 200", line) != 0) {
    gui_text_show(ERROR, "Not a valid Biograf 200 file.\n");
    return(1);
  if (g_ascii_strcasecmp("FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)", line) != 0)
  {
    gui_text_show(ERROR, "Unexpected atom format in Biograf file. Unable to continue.\n");
    return(2);
  }
*/

/* NEW - process (non-atom) header lines */
  do
    {
    if (fgetline(fp, line))
      {
      printf("Error: premature EOF\n");
      return(2);
      }
    }
  while (g_ascii_strncasecmp(line, "HETATM", 6) != 0);


/* all the next lines that start with HETATM are atom specifications */
  buff = tokenize(line, &num_tokens);

  while(g_ascii_strcasecmp("HETATM", *(buff)) == 0)
  {
    /* copy the first two characters */
    gchar *core_name = g_strndup(*(buff+2), 2);
    if (g_ascii_isdigit(core_name[1]))
    {
      g_free(core_name);
      core_name = g_strndup(*(buff+2), 1);
    }
    
    core = core_new(core_name, NULL, data);
    g_free(core_name);
    
    g_free(core->atom_type);
    core->atom_type = g_strdup(*(buff+9));
    
    /* default: all hydrogens can form hydrogen bonds, now with type information
      only H_HB types can  */
    if (g_ascii_strcasecmp("H_HB", core->atom_type) != 0)
    {  
      core->hydrogen_bond = FALSE;
    }
 
    /* Replace the 'normal' hydrogen type H_ by H_n for GULP/Dreiding */
    if (g_ascii_strcasecmp("H_", core->atom_type) == 0)
    {
      g_free(core->atom_type);
      core->atom_type = g_strdup("H_n");
    }
    
    g_free(core->atom_label);
    core->atom_label = NULL;
    /* lookup the atom type in the hash table */
    core->atom_label = g_strdup(g_hash_table_lookup(atom_types_and_labels, core->atom_type));
    /* if it's not found, than put the atom type (key) and the atom label (value) in the hash map */
    if (!core->atom_label)
    {
      g_hash_table_replace(atom_types_and_labels, g_strdup(core->atom_type), g_strdup(*(buff+2)));
      core->atom_label = g_strdup(*(buff+2));
    }
    
    /* use list appending, rather than prepending to get the right
       order in the list with respect to the atom_label */
    data->cores = g_slist_append(data->cores, core);
    
    core->x[0] = str_to_float(*(buff+6));
    core->x[1] = str_to_float(*(buff+7));
    core->x[2] = str_to_float(*(buff+8));
    
    /* the charge, based on Gasteiger is given as well */
    core->charge = str_to_float(*(buff+12));
    core->lookup_charge = FALSE;
    
    /* get the next line */
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
  }
  g_strfreev(buff);
 

/* FIXME - this isn't working - most likely because model_prep() hasn't been done yet */
#define DO_BGF_CONNECT 0
#if DO_BGF_CONNECT
  /* one empty line */
  fgetline(fp, line);
  
  /* next up is the bonding information
    it is formatted as FORMAT CONECT (a6,12i6)
    which means that up to twelve bonds can be accomodated for */
  /* CONECT record */
  buff = get_tokenized_line(fp, &num_tokens);
  /* ORDER record */
  buff2 = get_tokenized_line(fp, &num_tokens);
  
  while (g_ascii_strcasecmp(*(buff+0), "CONECT") == 0)
  {
    w = g_ascii_strtoull(*(buff+1), NULL, 10);
    core = g_slist_nth_data(data->cores, w - 1);
    for (i = 1; i < num_tokens - 1; ++i)
    {
      w = g_ascii_strtoull(*(buff + i + 1), NULL, 10);
      core_bond_to = g_slist_nth_data(data->cores, w - 1);
      w = g_ascii_strtoull(*(buff2 + i + 1), NULL, 10);
      w = ( (w == 1) ? BOND_SINGLE : ( (w == 2) ? BOND_DOUBLE : BOND_TRIPLE));
      connect_user_bond(core, core_bond_to, w, data);
    }
    /* free the strings before getting a new tokenized line */
    g_strfreev(buff);
    g_strfreev(buff2);
    buff = get_tokenized_line(fp, &num_tokens);
    buff2 = get_tokenized_line(fp, &num_tokens);
  }
  
  /* allocate 25 characters for each species line, that should be enough */
  gchar * species_string = g_new(gchar, g_hash_table_size(atom_types_and_labels) * 25 );
  /* get all species from the hash table */
  g_hash_table_foreach(atom_types_and_labels, species_line_from_hash_table, &species_string);
  
  if (data->gulp.species)
    g_free(data->gulp.species);
  data->gulp.species = species_string;
  
  /* all contained keys and values will be destroyed as well, because of the 
     g_hash_table_new_full init function */
  g_hash_table_destroy(atom_types_and_labels);
#endif

  
  /* got everything */
  data->num_asym = data->num_atoms = n;
  
  /* model setup */
  strcpy(data->filename, filename); 
  g_free(data->basename);
  data->basename = parse_strip(filename);
  
  model_prep(data);
  
  fclose(fp);
  return (0);
}

