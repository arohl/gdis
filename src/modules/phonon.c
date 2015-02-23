/*
Copyright (C) 2003 by Sean David Fleming

sean@power.curtin.edu.au

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
#include <math.h>

#include "../gdis.h"
#include "../coords.h"
#include "../matrix.h"
#include "../parse.h"
#include "../module.h"

/*****************************/
/* module configure function */
/*****************************/
gint module_libphonon_config(gpointer module)
{
module_label_set("phonon analysis", module);
module_symbol_add("show_phonon_components", module);
module_symbol_add("show_raman_active", module);
module_symbol_add("show_ir_active", module);
return(0);
}

/************************/
/* module functionality */
/************************/
void show_phonon_components(struct model_pak *model)
{
gdouble f, lc, lf, xc[3], xf[3];
gpointer ptr;
GSList *list;
struct core_pak *core;

/* checks */
if (!model)
  return;

/* get required mode to analyse */
ptr = g_slist_nth_data(model->phonons, model->current_phonon-1);
if (!ptr)
  return;
f = str_to_float(ptr);

printf("--------------------------------------------------------------------------\n");
printf("Mode: %d, Frequency = %f\n", model->current_phonon, f);
printf("--------------------------------------------------------------------------\n");
printf("  atom  |   len   |     x        y        z    |      a        b        c\n");
printf("--------------------------------------------------------------------------\n");

for (list=model->selection ; list; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;

/* get eigen-vector components */
  xc[0] = *((gdouble *) g_slist_nth_data(core->vibx_list, model->current_phonon-1));
  xc[1] = *((gdouble *) g_slist_nth_data(core->viby_list, model->current_phonon-1));
  xc[2] = *((gdouble *) g_slist_nth_data(core->vibz_list, model->current_phonon-1));
  ARR3SET(xf, xc);
  vecmat(model->ilatmat, xf);

  lc = VEC3MAG(xc);
  lf = VEC3MAG(xf);

  VEC3MUL(xc, 1.0/lc);
  VEC3MUL(xf, 1.0/lf);

  printf("%6s  | %7.4f | %7.2f  %7.2f  %7.2f  |  %7.2f  %7.2f  %7.2f\n",
          core->atom_label, lc, xc[0], xc[1], xc[2], xf[0], xf[1], xf[2]);

  }
}

/**********************************/
/* print non-zero intensity modes */
/**********************************/
void show_active(GSList *list, struct model_pak *model)
{
gint n;
gchar *freq;
gdouble intensity;
GSList *item;

n = 0;
for (item=list ; item ; item=g_slist_next(item))
  {
  intensity = str_to_float(item->data);

  if (intensity > 0.001)
    {
    freq = g_slist_nth_data(model->phonons, n);
    if (freq)
      printf(" %6d | %13s | %f\n", n+1, freq, intensity);
    }
  n++;
  } 
}

/***************************/
/* print non-zero IR modes */
/***************************/
void show_ir_active(struct model_pak *model)
{
if (!model)
  return;

printf("----------------------------------------\n");
printf(" mode # |   frequency   |  IR intensity\n");
printf("----------------------------------------\n");

show_active(model->ir_list, model);
}

/******************************/
/* print non-zero raman modes */
/******************************/
void show_raman_active(struct model_pak *model)
{
if (!model)
  return;

printf("------------------------------------------\n");
printf(" mode # |   frequency   | Raman intensity\n");
printf("------------------------------------------\n");
show_active(model->raman_list, model);
}
