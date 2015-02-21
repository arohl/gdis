/*
Copyright (C) 2004 by Craig Andrew James Fisher
Copyright (C) 2003 by Sean David Fleming

fisher@jfcc.or.jp

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
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "matrix.h"

#define PRECISION 1e-6 /* Precision of net charge */

/* globals for the model contents dialog */
GtkWidget *mcd_total, *mcd_charge, *mcd_molecules, *mcd_ions,
          *mcd_selected, *mcd_hidden;
/*****************************************/
/* updates the dialog for model contents */
/*****************************************/
#define DEBUG_ATOM_NOS 0
void atom_numbers_update(struct model_pak *data)
{
gint ntotal=0, nmols = 0, nions = 0, nhidden=0;
gchar *total, *charge, *molecules, *ions, *selected, *hidden;
GSList *list, *list2;
struct mol_pak *mol;
struct core_pak *core;
gdouble net_charge, molchg;

if (data)
  {
  net_charge = 0.0;
  for (list=data->moles ; list ; list=g_slist_next(list))
    {
    mol = (struct mol_pak *) list->data;
    ntotal += g_slist_length(mol->cores);

    molchg = 0.0;
    for (list2=mol->cores ; list2 ; list2=g_slist_next(list2))
      {
      core = (struct core_pak *) list2->data;
      molchg += atom_charge(core);
      if( core->status & HIDDEN )
        nhidden++;
      }
    if ( fabs(molchg) > 1e-6 )
      nions++; 
    else
      nmols++;
    net_charge += molchg;
    }
  if( fabs(net_charge) < PRECISION )
    net_charge = 0.0;

  selected = g_strdup_printf("%d",g_slist_length(data->selection));
  total = g_strdup_printf("%d", ntotal);
  charge = g_strdup_printf("%g", net_charge);
  ions = g_strdup_printf("%d", nions);
  molecules = g_strdup_printf("%d", nmols);
  hidden = g_strdup_printf("%d", nhidden);
  }
else
  {
  total = g_strdup("0");
  charge = g_strdup("0");
  ions = g_strdup("0");
  molecules = g_strdup("0");
  selected = g_strdup("0");
  hidden = g_strdup("0");
  }

/* entry updates */
gtk_entry_set_text(GTK_ENTRY(mcd_total), total);
gtk_entry_set_text(GTK_ENTRY(mcd_charge), charge);
gtk_entry_set_text(GTK_ENTRY(mcd_ions), ions);
gtk_entry_set_text(GTK_ENTRY(mcd_molecules), molecules);
gtk_entry_set_text(GTK_ENTRY(mcd_selected), selected);
gtk_entry_set_text(GTK_ENTRY(mcd_hidden), hidden);

gtk_entry_set_editable(GTK_ENTRY(mcd_total), FALSE);
gtk_entry_set_editable(GTK_ENTRY(mcd_charge), FALSE);
gtk_entry_set_editable(GTK_ENTRY(mcd_ions), FALSE);
gtk_entry_set_editable(GTK_ENTRY(mcd_molecules), FALSE);
gtk_entry_set_editable(GTK_ENTRY(mcd_selected), FALSE);
gtk_entry_set_editable(GTK_ENTRY(mcd_hidden), FALSE);

/* cleanup */
g_free(total);
g_free(molecules);
g_free(ions);
g_free(selected);
g_free(hidden);
}

/*********************************************/
/* Display total no and selected no of atoms */
/*********************************************/
void atom_numbers_box(GtkWidget *box)
{
GtkWidget *frame, *hbox, *vbox, *entry;

/* checks */
g_return_if_fail(box != NULL);

frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);

/* two column data display */
hbox = gtk_hbox_new(TRUE, 0);
gtk_container_add(GTK_CONTAINER(frame), hbox);

/* left vbox - titles */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_start(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);

/* TODO - put in a for loop? */

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Total atoms");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Net charge");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Ions");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Molecules");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Selected");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

entry = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry), "Hidden");
gtk_entry_set_editable(GTK_ENTRY(entry), FALSE);
gtk_box_pack_start(GTK_BOX(vbox), entry, FALSE, FALSE, 0);

/* right vbox - data */
vbox = gtk_vbox_new(TRUE, 0);
gtk_box_pack_end(GTK_BOX(hbox), vbox, TRUE, TRUE, 0);

mcd_total = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_total, FALSE, FALSE, 0);

mcd_charge = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_charge, FALSE, FALSE, 0);

mcd_ions = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_ions, FALSE, FALSE, 0);

mcd_molecules = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_molecules, FALSE, FALSE, 0);

mcd_selected = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_selected, FALSE, FALSE, 0);

mcd_hidden = gtk_entry_new();
gtk_box_pack_start(GTK_BOX(vbox), mcd_hidden, FALSE, FALSE, 0);

/* set boxes */
atom_numbers_update(NULL);
}
