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

#include "gdis.h"
#include "coords.h"
#include "interface.h"

extern struct sysenv_pak sysenv;

struct undo_pak
{
gint (*function) (gpointer, gpointer);
GSList *pointers;
};

/********************************/
/* initialize model's undo list */
/********************************/
void undo_init(struct model_pak *model)
{
model->undo_list = NULL;
}

/*************************************************/
/* removes all registered undo items for a model */
/*************************************************/
void undo_free(struct model_pak *model)
{
GSList *list;
struct undo_pak *undo;

g_assert(model != NULL);

list = model->undo_list;
while (list)
  {
  undo = list->data;
  list = g_slist_next(list);

  free_slist(undo->pointers);
  g_free(undo);
  }

g_slist_free(model->undo_list);
}

/***********************************/
/* record an undo item for a model */
/***********************************/
void undo_register(struct model_pak *model, gpointer f, GSList *p)
{
struct undo_pak *undo;

g_assert(model != NULL);

undo = g_malloc(sizeof(struct undo_pak));
undo->function = f;
undo->pointers = p;

model->undo_list = g_slist_prepend(model->undo_list, undo);
}

/********************************************/
/* invoke the topmost undo item for a model */
/********************************************/
void undo_single(struct model_pak *model)
{
GSList *list;
struct undo_pak *undo;

if (!model)
  return;

list = model->undo_list;
if (list)
  {
  undo = list->data;
  g_assert(undo != NULL);

  if (undo->function)
    undo->function(model, undo->pointers);

  model->undo_list = g_slist_remove(model->undo_list, undo);
  }
}

/*************************************************/
/* undo most recent for the current active model */
/*************************************************/
void undo_active(void)
{
undo_single(sysenv.active_model);
gui_refresh(GUI_CANVAS);
}

