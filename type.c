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
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "type.h"
#include "type_pak.h"
#include "task.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* typing ruleset construction */
/* FIXME - associate with model? globally? either/or? eg like element data */

/****************************************/
/* print the contents of type structure */
/****************************************/
void type_print(gpointer data)
{
struct type_pak *type = data;

printf(" *** Typing rule: %p\n", type);
printf("  FF assignment = %s\n", type->ff);
printf("  Rule list, size = %d\n", g_slist_length(type->rules));

}

/***********************************/
/* create a new atom typing object */
/***********************************/
gpointer type_new(void)
{
struct type_pak *type;

type = g_malloc(sizeof(struct type_pak));

type->set_ff = FALSE;
type->ff = NULL;

type->set_charge = FALSE;
type->charge = 0.0;

type->rules = NULL;

return(type);
}

/************************/
/* free a typing object */
/************************/
void type_free(gpointer data)
{
struct type_pak *type = data;

g_assert(type != NULL);

g_free(type->ff);

free_slist(type->rules);

g_free(type);
}

/*************************************************/
/* search for undefined types in a list of cores */
/*************************************************/
gint type_check_list(GSList *cores)
{
gint n=0;
struct core_pak *core;
GSList *list;

for (list=cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
/* NULL or '?' */
  if (!core->atom_type)
    n++;
  }
return(n);
}

/***********************************/
/* configure forcefield assignment */
/***********************************/
void type_ff_set(gint flag, const gchar *ff, gpointer data)
{
struct type_pak *type = data;

g_assert(type != NULL);

type->set_ff = flag;
type->ff = g_strdup(ff);
}

/*******************************/
/* configure charge assignment */
/*******************************/
void type_charge_set(gint flag, gdouble charge, gpointer data)
{
struct type_pak *type = data;

g_assert(type != NULL);

type->set_charge = flag;
type->charge = charge;
}

/***************************************/
/* add a rule to an atom typing object */
/***************************************/
/* TODO - more flexible eg not only elem, but label matching? */
void type_rule_add(gint level, gint count, const gchar *elem, gpointer data)
{
struct type_pak *type = data;
struct rule_pak *rule;

rule = g_malloc(sizeof(struct rule_pak));

rule->atom_number = elem_symbol_test(elem);
rule->count = count;
rule->level = level;

type->rules = g_slist_prepend(type->rules, rule);
}

/*****************************************/
/* test a rule against a particular atom */
/*****************************************/
gint type_rule_match(struct rule_pak *rule, struct core_pak *core)
{
gint count;
GSList *list, *neighbours;
struct core_pak *core2;

switch (rule->level)
  {
  case 0:
    if (rule->atom_number == core->atom_code)
      return(TRUE);
    else
      return(FALSE);

  case 1:
    neighbours = connect_neighbours(core);
    count = 0;
    for (list=neighbours ; list ; list=g_slist_next(list))
      {
      core2 = list->data;
      if (rule->atom_number == core2->atom_code)
        count++;
      }
    if (count >= rule->count)
      return(TRUE);
    else
      return(FALSE);

  default:
    printf("FIXME - rule matching at level %d not implemented yet.\n", rule->level);
  }
return(FALSE);
}

/******************************************/
/* apply a typing rule to a list of atoms */
/******************************************/
gint type_apply(gpointer data, GSList *cores)
{
gint match;
GSList *list1, *list2;
struct core_pak *core;
struct type_pak *type = data;

g_assert(type != NULL);

/* foreach atom - check it satisfies all rules in the type's rule list */
for (list1=cores ; list1 ; list1=g_slist_next(list1))
  {
  core = list1->data;

/* rule match */
  match = TRUE;
  for (list2=type->rules ; list2 ; list2=g_slist_next(list2))
    {
/* check for failed rule (only need 1 failure to terminate) */
    if (!type_rule_match(list2->data, core))
      {
      match = FALSE; 
      break;
      }
    }

  if (match)
    {
/*
printf("core %p : assigning type %s\n", core, type->ff);
*/
    g_free(core->atom_type);
    core->atom_type = g_strdup(type->ff);
    }
  }

/* TODO - return the number matched? */
return(0);
}

/**************************/
/* remove all atom typing */
/**************************/
void type_clear(struct model_pak *model)
{
GSList *clist;
struct core_pak *core;

g_assert(model != NULL);

for (clist=model->cores ; clist ; clist=g_slist_next(clist))
  {
  core = clist->data;
  if (core->atom_type)
    {
    g_free(core->atom_type);
    core->atom_type = NULL;
    }
  }
}

/*************************************************************************/
/* get the number of bonded atoms of type atom_code attached to the core */
/*************************************************************************/
#define DEBUG_TYPE_CORE_HAS_BOND 0
gint type_core_has_bond(gint atom_code, struct core_pak *core) 
{
gint match=0;
GSList *blist;
struct bond_pak *bond;
struct core_pak *core1, *core2;

#if DEBUG_TYPE_CORE_HAS_BOND
printf("bond search [%s] - [%s]:\n", core->label, elements[atom_code].symbol);
#endif

for (blist=core->bonds ; blist ; blist=g_slist_next(blist))
  {
  bond = blist->data;

  if (bond->type == BOND_HBOND)
    continue;

  core1 = bond->atom1;
  core2 = bond->atom2;

  if (core1 == core)
    {
    if (core2->atom_code == atom_code)
      match++; 
    }
  else
    {
    if (core1->atom_code == atom_code)
      match++; 
    }
  }

#if DEBUG_TYPE_CORE_HAS_BOND
printf("%d matches found.\n", match);
#endif

return(match);
}

/********************************************/
/* assign CVFF atom labels to a given model */
/********************************************/
void type_as_cvff(struct model_pak *model)
{
gint nb, nc;
GSList *clist, *mlist;
struct core_pak *core;
struct mol_pak *mol;

g_assert(model != NULL);

/* wipe any pre-existing typing */
type_clear(model);

/* loop over molecules -> can get some useful info (eg charge/size) */
for (mlist=model->moles ; mlist ; mlist=g_slist_next(mlist))
  {
  mol = mlist->data;
  nc = g_slist_length(mol->cores);

/* TODO - get molecule charge */

  for (clist=mol->cores ; clist ; clist=g_slist_next(clist))
    {
    core = clist->data;
    nb = g_slist_length(core->bonds);

/* skip any that have already been assigned */
/* TODO - convenient if (eg) we latch onto an o with 2 h's (mol size == 3) -> water */
    if (core->atom_type)
      continue;

    switch(core->atom_code)
      {
/* hydrogen */
      case 1:
        if (type_core_has_bond(8, core))
          {
/* FIXME - not a very good test for water... */
          if (nc == 3)
            core->atom_type = g_strdup("h*");
          else 
            core->atom_type = g_strdup("ho");
          }
        else
          core->atom_type = g_strdup("h");
        break;

/* carbon */
      case 6:
        switch(nb)
          {
          case 4:
            core->atom_type = g_strdup("cg");
            break; 
 
          case 3:
/* check -c-   
          |    
          h  */
/* carboxylate */
            if (type_core_has_bond(8, core) == 2)
              core->atom_type = g_strdup("c-");
            else
              core->atom_type = g_strdup("cp");
            break;

          case 2:
/* check -c-o , o-c-o */
            core->atom_type = g_strdup("c");
            break;

          case 1:
/* always a triple? (eg like co) */
            core->atom_type = g_strdup("ct");
            break;

          default:
            core->atom_type = g_strdup("c");
          }
        break;

/* oxygen */
      case 8:
/* TODO - better scheme needed - if bonded to carbon - need to check that carbon */
/* TODO - initially set all types to NULL & then while(!type) ... */
        switch (nb)
          {
          case 2:
            if (type_core_has_bond(1, core) == 2)
              core->atom_type = g_strdup("o*");
            else
              {
              if (type_core_has_bond(1, core) == 1)
                core->atom_type = g_strdup("oh");
              else
                core->atom_type = g_strdup("o");
              }
            break;

          case 1:
            if (type_core_has_bond(1, core))
              core->atom_type = g_strdup("oh");
            else
              core->atom_type = g_strdup("o");
            break;

          default:
            core->atom_type = g_strdup("o");
            break;
          }
        break;

      default:
/* FIXME - is just the element symbol sufficient? */
        core->atom_type = g_strdup_printf("%s", elements[core->atom_code].symbol);
        break;
      }
    }
  }
}

/**********************************************/
/* duplicate a model's cores & shells exactly */
/**********************************************/
/* TODO - put in coords.c ? */
void copy_atom_data(struct model_pak *src, struct model_pak *dest)
{
GSList *list;
struct core_pak *core;

/* copy periodicity info */
dest->periodic = src->periodic;
dest->fractional = src->fractional;
memcpy(dest->pbc, src->pbc, 6*sizeof(gdouble));

/* copy cores & shells */
for (list=src->cores ; list ; list=g_slist_next(list))
  {
  core = dup_core(list->data);
  dest->cores = g_slist_prepend(dest->cores, core);

  if (core->shell)
    dest->shels = g_slist_prepend(dest->shels, core->shell);
  }

/* ensure cores & shells are in the same order */
dest->cores = g_slist_reverse(dest->cores);
dest->shels = g_slist_reverse(dest->shels);
}

/**********************************************/
/* assign atom charges (EEM) to a given model */
/**********************************************/
void type_eem(struct model_pak *model)
{
/* TODO */
}

/**********************************************/
/* assign atom charges (QEq) to a given model */
/**********************************************/
void type_qeq(struct model_pak *model)
{
gint num_tokens;
gchar line[LINELEN];
gchar *inp, *out, **buff;
GSList *list;
struct core_pak *core, *parent;
struct model_pak temp;
FILE *fp;

/* checks */
g_assert(model != NULL);

printf("Calculating charges using QEq...\n");

/* temp model for QEq calc */
model_init(&temp);
copy_atom_data(model, &temp);
/* setup for GULP QEq calc */
temp.gulp.run = E_SINGLE;
temp.gulp.qeq = TRUE;
temp.gulp.method = model->gulp.method;
/* coulomb must be disabled for QEq */
temp.gulp.coulomb = NOBUILD;
/* use symmetry for speed */
temp.sginfo.spacename = g_strdup(model->sginfo.spacename);
temp.sginfo.spacenum = model->sginfo.spacenum;

/* get temporary input/output filenames */
inp = gun("gin");
out = gun("got");

/* use GULP for the QEq calc */
write_gulp(inp, &temp);
exec_gulp(inp, out);

/* assign charges */
fp = fopen(out, "rt");
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "  Final charges from QEq", 24) == 0)
    {
/* proceed until 1st atom */
    buff = get_tokenized_line(fp, &num_tokens);
    while (buff)
      {
      if (g_ascii_isdigit(**(buff+0)) && g_ascii_isdigit(**(buff+1)))
        break;
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);
      }

/* got 1st valid atom */
    list = model->cores;
    while (buff)
      {
/* check */
      if (num_tokens < 3)
        break;

/* current core to update */
      while (list)
        {
        core = list->data;

/* NB: we're using the asymetric unit (if exists) */
        if (core->primary)
          break;
        else
          list = g_slist_next(list);
        }
      if (list)
        list = g_slist_next(list);
      else
        break;

/*
printf("%s : (%s) %s\n", core->label, *(buff+1), *(buff+2));
*/

/* TODO - what about shells??? */

      core->charge = str_to_float(*(buff+2));
      core->lookup_charge = FALSE;

      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);
      }
    }
  }
fclose(fp);

/* apply charges to symmetry related atoms as well */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->primary)
    continue;

  g_assert(core->primary_core != NULL);
  parent = core->primary_core;

  core->charge = parent->charge;
  core->lookup_charge = FALSE;
  }
calc_emp(model);

/* cleanup (TODO - delete inp/out) */
model_free(&temp);
}

#define STRICT_DREIDING_TYPE_CHECK 0

void type_dreiding_gasteiger(struct model_pak *model, gint mode)
{
  gchar *temp_in_file, *temp_out_file;
  gchar *cmd;
  
  GSList *atom_list, *typed_atom_list;
  
  struct core_pak *core, *typed_core;
  struct model_pak temp_model;
  
  /* be sure that babel is installed */
  if (!sysenv.babel_path)
  {
    gui_text_show(ERROR, "The babel program was not found in your path.\n");
    return;
  }
  
  /* define temporary files */
/*
#ifdef __WIN32
  temp_in_file = "C:\\babel_in_temp.xyz";
  temp_out_file = "C:\\babel_out_temp.bgf";
#else
  temp_in_file = "/tmp/babel_in_temp.xyz";
  temp_out_file = "/tmp/babel_out_temp.bgf";
#endif
  write_xyz(temp_in_file, model);
  cmd = g_strdup_printf("%s -ixyz %s -obgf %s", sysenv.babel_path, temp_in_file, temp_out_file);
*/


/* NEW - using PDB as temp file, since it has connect data */
/* which babel may not generate correctly (at all?) for XYZ */
#ifdef __WIN32
  temp_in_file = "C:\\babel_in_temp.pdb";
  temp_out_file = "C:\\babel_out_temp.bgf";
#else
  temp_in_file = "/tmp/babel_in_temp.pdb";
  temp_out_file = "/tmp/babel_out_temp.bgf";
#endif

  write_pdb(temp_in_file, model);
  cmd = g_strdup_printf("%s -ipdb %s -obgf %s", sysenv.babel_path, temp_in_file, temp_out_file);

  
  system(cmd);
  g_free(cmd);
  
  /* create a new temporary model */
  model_init(&temp_model);
  
  /* read the translated result */
  read_bgf(temp_out_file, &temp_model);  
  
#ifdef STRICT_DREIDING_TYPE_CHECK
  g_return_if_fail(g_slist_length(model->cores) == g_slist_length(temp_model.cores));
#endif
  /* dual list traversal */
  for (atom_list = model->cores, typed_atom_list = temp_model.cores;
       atom_list; 
       atom_list = atom_list->next, typed_atom_list = typed_atom_list->next)
  {
    core = atom_list->data;
    typed_core = typed_atom_list->data;
#if STRICT_DREIDING_TYPE_CHECK
    g_return_if_fail(core->atom_code == typed_core->atom_code);
    /* can't check location or label, as the location can be fractional and 
      the label will almost certainly be different. */
#endif
    if (mode == DREIDING)
    {
      g_free(core->atom_type); 
      core->atom_type = g_strdup(typed_core->atom_type);

/*
      g_free(core->atom_label);
      core->atom_label = g_strdup(typed_core->atom_label);
*/

      if (model->gulp.species)
        g_free(model->gulp.species);
      model->gulp.species = g_strdup(temp_model.gulp.species);
    }
    if (mode == GASTEIGER) 
    {
      core->charge = typed_core->charge;
      core->lookup_charge = FALSE;
    }
  }  
  model_free(&temp_model);
  
  unlink(temp_in_file);
  unlink(temp_out_file);
}

/*********************************/
/* type model according to label */
/*********************************/
void type_model(const gchar *type, struct model_pak *model)
{
if (g_ascii_strncasecmp("CVFF", type, 4) == 0)
  type_as_cvff(model);
if (g_ascii_strncasecmp("QEq", type, 3) == 0)
  type_qeq(model);
if (g_ascii_strncasecmp("Drei", type, 4) == 0)
  type_dreiding_gasteiger(model, DREIDING);
if (g_ascii_strncasecmp("Gast", type, 4) == 0)
  type_dreiding_gasteiger(model, GASTEIGER);
}
