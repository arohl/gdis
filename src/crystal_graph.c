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
#include <stdlib.h>
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
#include "type.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

GSList * create_pair(struct model_pak * model, struct mol_pak * mol1, struct mol_pak * mol2, gint x, gint y, gint z)
{
  GSList * core_iter = NULL, *mol2_copy = NULL;
  struct core_pak * core;
  if (mol1 == mol2 && 
      x == 0 &&
      y == 0 &&
      z == 0)
    return NULL;
  
  /* copy cores */
  mol2_copy = dup_core_list(mol2->cores);
  
  /* create a symmetry-related molecule and convert its fractional coordinates to cartesian*/
  for (core_iter = mol2_copy; core_iter; core_iter = core_iter->next)
  {
    core = core_iter->data;
    core->primary = TRUE; //have to make primary, otherwise it will not be printed in gulp input file!
    core->x[0] += x;
    core->x[1] += y;
    core->x[2] += z;
    vecmat(model->latmat, core->x);
  }
  
  /* add copies of the mol1 cores and convert their fractional coordinates as well */
  for (core_iter = mol1->cores; core_iter; core_iter = core_iter->next)
  {
    core = dup_core(core_iter->data);
    core->primary = TRUE;
    vecmat(model->latmat, core->x);
    mol2_copy = g_slist_prepend(mol2_copy, core);
  }
  /* return the result */
  return mol2_copy;
}

/* some tree functions for setting up a tree with energies as keys */

GMemChunk *key_chunk;
//GTree *tree;

gint energy_comp(gconstpointer a_ptr, gconstpointer b_ptr, gpointer ignored)
{
  gdouble a, b;
  a = *(gdouble *)a_ptr;
  b = *(gdouble *)b_ptr;
  if (a < b)
    return 1;
//  if (fabs(a - b) < 1e-10) //no equal keys, every thing must end up in the tree
//    return 0;
  return -1; /* a >= b */ 
}

void free_key(gpointer key)
{
  g_mem_chunk_free(key_chunk, key);
}

void free_value(gpointer value)
{
  g_free((gchar *)value);
}

gboolean print_node(gpointer key, gpointer value, gpointer ignored)
{
  g_print("%-20s have interaction energy of %3.5f kJ/mol\n", (gchar *)value, *(gdouble *)key);
  return FALSE;
}

gdouble run_gulp(GSList *cores, gchar *species)
{
  struct model_pak *new_model, *result_model;  
  new_model = g_malloc(sizeof(struct model_pak));
  result_model = g_malloc(sizeof(struct model_pak));
    
  model_init(new_model);
  model_init(result_model);
  new_model->cores = cores;
  /* init gulp parameters */
  new_model->gulp.run = E_SINGLE;
  new_model->gulp.method = CONP;
  new_model->gulp.species = g_strdup(species);
  new_model->gulp.libfile = g_strdup("dreiding.lib");
  new_model->gulp.no_exec = TRUE; //for now we only write the input file
  new_model->gulp.temp_file = g_strdup("cgf_pair.gin");
  new_model->gulp.out_file = g_strdup("cgf_pair.got");
  new_model->gulp.coulomb = MOLE;
  new_model->periodic = 0;
  
  write_gulp(new_model->gulp.temp_file, new_model);
  gchar * cmd = g_strdup_printf("%s < %s > %s", sysenv.gulp_path, 
                                new_model->gulp.temp_file, new_model->gulp.out_file);
  system(cmd);
  g_free(cmd);
  
  read_gulp_output(new_model->gulp.out_file, result_model);
  gdouble result = result_model->gulp.energy;
  model_free(new_model);
  model_free(result_model);

  return result;
}

gint calculate_crystal_graph(struct model_pak * model)
{
  gint m1 = 0, m2 = 0;

  GTree *tree;

  struct mol_pak *mol1, *mol2;
  struct core_pak *core;
  GSList *pair = NULL;
  GSList *mol_iter= NULL, *mol_iter2 = NULL, *mol_copy = NULL, *core_iter = NULL;
  gint x = 0, y = 0, z = 0;
  
  gint a_images =(gint) model->monty.image_x + 0.5;
  gint b_images =(gint) model->monty.image_y + 0.5;
  gint c_images =(gint) model->monty.image_z + 0.5;
  
  /* checks */ 
  g_return_val_if_fail(sysenv.gulp_path != NULL, 4);
  g_return_val_if_fail(model != NULL, 1);
  g_return_val_if_fail(a_images > 0 &&  b_images > 0 && c_images > 0, 2);
  if (model->periodic != 3)
  {
    gui_text_show(ERROR, "Can not calculate a crystal graph for non crystalline models");
    return 3;
  }
  gint nr_molecules = g_slist_length(model->moles);
  /* initialize the key chunk to hold the maximum number of  */
  key_chunk = g_mem_chunk_create(gint, 
                                 nr_molecules * nr_molecules *(2*a_images+1)*(2*b_images+1)*(2*c_images+1),
                                 G_ALLOC_AND_FREE);
  
  tree = g_tree_new_full(energy_comp, 
                         NULL,
                         free_key,
                         free_value);  
  
  /* type to get the species line */
  type_dreiding_gasteiger(model, DREIDING);
  
  /* calculate the single energies */
  gdouble * single_energies = g_new0(gdouble, nr_molecules);
  
  for (mol_iter = model->moles; mol_iter; mol_iter = mol_iter->next)
  {
    mol1 = mol_iter->data;
    mol_copy = dup_core_list(mol1->cores);
    /* add copies of the mol1 cores and convert their fractional coordinates as well */
    for (core_iter = mol_copy; core_iter; core_iter = core_iter->next)
    {
      core = core_iter->data;
      core->primary = TRUE;
      vecmat(model->latmat, core->x);
    }
    single_energies[m1] = run_gulp(mol_copy, model->gulp.species);
    single_energies[m1++] *= 96.48474;
  }

  m1 = 0;
  /* iterate over all molecule combinations */
  for (mol_iter = model->moles; mol_iter; mol_iter = mol_iter->next)
  {
    mol1 = mol_iter->data;
    m2 = 0;
    for (mol_iter2 = model->moles; mol_iter2; mol_iter2 = mol_iter2->next)
    {
      mol2 = mol_iter2->data;
      /* create all x, y and z images
         this can be done smarter, as all bonds are now calculated twice */
      for (x = -a_images; x <= a_images; ++x)
      {
        for (y = -b_images; y <= b_images; ++y)
        {
          for (z = -c_images; z <= c_images; ++z)
          {
            pair = create_pair(model, mol1, mol2, x, y, z);
            if (pair != NULL)
            {              
              gdouble *key_ptr;
              key_ptr = g_chunk_new(gdouble, key_chunk);
              *key_ptr = run_gulp(pair, model->gulp.species);
              *key_ptr *= 96.48474; //eV --> kJ / mol
              *key_ptr -= single_energies[m1];
              *key_ptr -= single_energies[m2];
              gchar * combination = g_strdup_printf(" %i - %i [ %i %i %i ]", m1+1, m2+1, x,y,z);
              g_tree_insert(tree, key_ptr, combination);
            }
          }
        }
      }
      m2++;
    }
    m1++;
  }
  /* new iterate the results from the tree */
  g_tree_foreach(tree, print_node, NULL);
  g_tree_destroy(tree);
  return(0);
}


