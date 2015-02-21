/*
Copyright (C) 2000 by Sean David Fleming

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
#include <unistd.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "parse.h"
#include "matrix.h"
#include "zmatrix.h"
#include "zmatrix_pak.h"
#include "gui_shorts.h"
#include "dialog.h"
#include "interface.h"
#include "zone.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

GtkWidget *zmat_combo, *zmat_entry;

/******************************************************/
/* process text window in order to update the zmatrix */
/******************************************************/
void gui_zmatrix_text_parse(GtkWidget *w, gpointer dialog)
{
gint n, tokens;
gchar *text, **buff;
GSList *list;
GtkWidget *view;
GtkTextBuffer *buffer;
GtkTextIter start, end;
struct model_pak *model;
struct zmat_pak *zmat;
struct zval_pak *zval;

/* acquire model data */
model = dialog_model(dialog);
g_assert(model != NULL);

zmat = model->zmatrix;
g_assert(zmat != NULL);

/* acquire dialog data */
view = dialog_child_get(dialog, "view1");
g_assert(GTK_IS_TEXT_VIEW(view));

buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
g_assert(GTK_IS_TEXT_BUFFER(buffer));

n=0;
for (list=zmat->zlines ; list ; list=g_slist_next(list))
  {
  zval = list->data;

/* extract desired line */
  gtk_text_buffer_get_iter_at_line_offset(buffer, &start, n, 0);
  gtk_text_buffer_get_iter_at_line_offset(buffer, &end, n+1, 0);
  text = gtk_text_buffer_get_slice(buffer, &start, &end, TRUE);

/* separate data items */
  buff = tokenize(text, &tokens);

/* set text entry values */
  if (tokens > 7)
    {
    if (str_is_float(*(buff+5)))
      zval->value[0] = str_to_float(*(buff+5));

    if (str_is_float(*(buff+6)))
      zval->value[1] = str_to_float(*(buff+6));

    if (str_is_float(*(buff+7)))
      zval->value[2] = str_to_float(*(buff+7));

/*
printf("+[%d]", n+1);
P3VEC(" ", zval->value);
*/

    }

  g_strfreev(buff);

  n++;
  }

}

/******************************************/
/* redo the zmatrix molecule from scratch */
/******************************************/
void gui_zmatrix_model_update(GtkWidget *w, gpointer dialog)
{
GSList *list;
struct zmat_pak *zmat;
struct model_pak *model;

/* acquire model data */
model = dialog_model(dialog);
g_assert(model != NULL);


/* CURRENT - acquire any changes to zmatrix text */
/*
gui_zmatrix_text_parse(NULL, dialog);
*/


zmat = model->zmatrix;

/* delete old zmatrix cores */
for (list=zmat->zcores ; list ; list=g_slist_next(list))
  delete_core(list->data);
delete_commit(model);
g_slist_free(zmat->zcores);
zmat->zcores = NULL;

/* recompute zmatrix core list */
zmat_process(model->zmatrix, model);

/* zmatrix cores always created as cartesians - but if the */
/* model is periodic - they are expected to be fractional */
if (!model->fractional && model->periodic)
  {
  for (list=zmat->zcores ; list ; list=g_slist_next(list))
    {
    struct core_pak *core = list->data;
    vecmat(model->ilatmat, core->x);
    }
  }

/* assign element data to cores */
for (list=zmat->zcores ; list ; list=g_slist_next(list))
  elem_init(list->data, model);

/* refresh coords/connectivity */
zone_init(model);
coords_compute(model);
connect_refresh(model);

gui_refresh(GUI_CANVAS);
}

/*****************************************************/
/* update zmatrix hash table with gui variable value */
/*****************************************************/
void gui_zmatrix_value_update(GtkWidget *w, gpointer data)
{
const gchar *key, *value;
struct zmat_pak *zmat = data;

key = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(zmat_combo)->entry));
value = gtk_entry_get_text(GTK_ENTRY(zmat_entry));

g_hash_table_insert(zmat->vars, g_strdup(key), g_strdup(value));
}

/********************************************************/
/* change current variable value to match variable name */
/********************************************************/
void gui_zmatrix_value_change(GtkWidget *w, gpointer data)
{
const gchar *text, *value;
struct zmat_pak *zmat = data;

text = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(zmat_combo)->entry));

value = g_hash_table_lookup(zmat->vars, text);

if (value)
  gtk_entry_set_text(GTK_ENTRY(zmat_entry), value);
}

/****************************************/
/* write all variable names into a list */
/****************************************/
void gui_zmatrix_vars(gpointer key, gpointer val, gpointer data)
{
GList **list = data;

*list = g_list_append(*list, key);
}

/*****************************/
/* update total zmatrix text */
/*****************************/
void gui_zmatrix_text_update(GtkWidget *w, gpointer dialog)
{
gint i, n;
GSList *list;
GString *buff;
GtkWidget *view;
struct model_pak *model;
struct zmat_pak *zmat;
struct zval_pak *zval;

model = dialog_model(dialog);
g_assert(model != NULL);

zmat = model->zmatrix;
g_assert(zmat != NULL);

view = dialog_child_get(dialog, "view1");
g_assert(view != NULL);

/* generate the zmatrix lines for display */
n=1;
buff = g_string_new(NULL);
for (list=zmat->zlines ; list ; list=g_slist_next(list))
  {
  zval = list->data;

  g_string_append_printf(buff, "[%d]  %s  %d %d %d", n++, zval->elem,
                         zval->connect[0], zval->connect[1], zval->connect[2]);
  for (i=0 ; i<3 ; i++)
    {
    if (zval->name[i])
      g_string_append_printf(buff, "  %-9s", (gchar *) zval->name[i]);
    else
      g_string_append_printf(buff, "  %-9.4f", zval->value[i]);
    }
  g_string_append_printf(buff, "\n");
  }

gtk_text_buffer_set_text(gtk_text_view_get_buffer(GTK_TEXT_VIEW(view)), buff->str, -1);

g_string_free(buff, TRUE);
}

/************************************/
/* modify lines in the zmatrix text */
/************************************/
void gui_zmatrix_line_changed(GtkWidget *w, gpointer dialog)
{
gint n;
gchar *line;
const gchar *text1, *text2, *text3;
gpointer entry1, entry2, entry3;
GtkWidget *view, *spin;
GtkTextBuffer *buffer;
GtkTextIter start, end;
struct model_pak *model;
struct zmat_pak *zmat;
struct zval_pak *zval;

/* acquire widget information */
view = dialog_child_get(dialog, "view1");
g_assert(GTK_IS_TEXT_VIEW(view));

spin = dialog_child_get(dialog, "spin1");
g_assert(GTK_IS_SPIN_BUTTON(spin));

buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
g_assert(GTK_IS_TEXT_BUFFER(buffer));

n = SPIN_IVAL(GTK_SPIN_BUTTON(spin));

entry1 = dialog_child_get(dialog, "entry1");
entry2 = dialog_child_get(dialog, "entry2");
entry3 = dialog_child_get(dialog, "entry3");

g_assert(GTK_IS_ENTRY(entry1));
g_assert(GTK_IS_ENTRY(entry2));
g_assert(GTK_IS_ENTRY(entry3));

text1 = gtk_entry_get_text(GTK_ENTRY(entry1));
text2 = gtk_entry_get_text(GTK_ENTRY(entry2));
text3 = gtk_entry_get_text(GTK_ENTRY(entry3));


/* acquire model information */
model = dialog_model(dialog);
g_assert(model != NULL);

zmat = model->zmatrix;
g_assert(zmat != NULL);

/* form a new zmatrix text line */
zval = g_slist_nth_data(zmat->zlines, n-1);

/* out of range check */
if (!zval)
  return;

line = g_strdup_printf("[%d]  %s  %d %d %d  %-9s  %-9s  %-9s\n", n, zval->elem,
                  zval->connect[0], zval->connect[1], zval->connect[2],
                  text1, text2, text3);

/* start of line */
gtk_text_buffer_get_iter_at_line_offset(buffer, &start, n-1, 0);
/* effective end of line */
/* FIXME - actually the start of the next line - will this cause 1 char offset problems? */
gtk_text_buffer_get_iter_at_line_offset(buffer, &end, n, 0);
gtk_text_buffer_delete(buffer, &start, &end);


/* update appropriate zmatrix value with entry value */
if (entry1 == w)
  {
  if (str_is_float(text1))
    zval->value[0] = str_to_float(text1);
  }
if (entry2 == w)
  {
  if (str_is_float(text2))
    zval->value[1] = str_to_float(text2);
  }
if (entry3 == w)
  {
  if (str_is_float(text3))
    zval->value[2] = str_to_float(text3);
  }


gtk_text_buffer_insert_with_tags_by_name(buffer, &start, line, strlen(line), "bg_changed", NULL);

/*
if (changed)
  gtk_text_buffer_insert_with_tags_by_name(buffer, &start, line, strlen(line), "bg_changed", NULL);
else
  gtk_text_buffer_insert(buffer, &start, line, strlen(line));
*/

}

/**************************************/
/* aqcuire a zmatrix line for editing */
/**************************************/
void gui_zmat_line_update(GtkWidget *w, gpointer dialog)
{
gint n, tokens;
gchar *text1, **buff;
gpointer entry1, entry2, entry3;
GtkWidget *view, *spin;
GtkTextBuffer *buffer;
GtkTextIter start, end;

/* acquire dialog data */
view = dialog_child_get(dialog, "view1");
g_assert(GTK_IS_TEXT_VIEW(view));

buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
g_assert(GTK_IS_TEXT_BUFFER(buffer));

spin = dialog_child_get(dialog, "spin1");
g_assert(GTK_IS_SPIN_BUTTON(spin));

n = SPIN_IVAL(GTK_SPIN_BUTTON(spin));

entry1 = dialog_child_get(dialog, "entry1");
entry2 = dialog_child_get(dialog, "entry2");
entry3 = dialog_child_get(dialog, "entry3");

/* extract desired line */
gtk_text_buffer_get_iter_at_line_offset(buffer, &start, n-1, 0);
gtk_text_buffer_get_iter_at_line_offset(buffer, &end, n, 0);
text1 = gtk_text_buffer_get_slice(buffer, &start, &end, TRUE);

/* separate data items */
buff = tokenize(text1, &tokens);

/* stop entry changes from generating an entry changed */
/* signal since we're forcibly setting the entry */
g_signal_handlers_block_by_func(entry1, gui_zmatrix_line_changed, dialog);
g_signal_handlers_block_by_func(entry2, gui_zmatrix_line_changed, dialog);
g_signal_handlers_block_by_func(entry3, gui_zmatrix_line_changed, dialog);

/* set text entry values */
if (tokens > 7)
  {
  gtk_entry_set_text(GTK_ENTRY(entry1), *(buff+5));
  gtk_entry_set_text(GTK_ENTRY(entry2), *(buff+6));
  gtk_entry_set_text(GTK_ENTRY(entry3), *(buff+7));

/* disallow editing of variable names */
  if (str_is_float(*(buff+5)))
    gtk_widget_set_sensitive(GTK_WIDGET(entry1), TRUE);
  else
    gtk_widget_set_sensitive(GTK_WIDGET(entry1), FALSE);

  if (str_is_float(*(buff+6)))
    gtk_widget_set_sensitive(GTK_WIDGET(entry2), TRUE);
  else
    gtk_widget_set_sensitive(GTK_WIDGET(entry2), FALSE);

  if (str_is_float(*(buff+7)))
    gtk_widget_set_sensitive(GTK_WIDGET(entry3), TRUE);
  else
    gtk_widget_set_sensitive(GTK_WIDGET(entry3), FALSE);
  }
else
  {
  gtk_entry_set_text(GTK_ENTRY(entry1), "");
  gtk_entry_set_text(GTK_ENTRY(entry2), "");
  gtk_entry_set_text(GTK_ENTRY(entry3), "");
  }

g_strfreev(buff);

/* allow normal entry change processing to continue */
g_signal_handlers_unblock_by_func(entry1, gui_zmatrix_line_changed, dialog);
g_signal_handlers_unblock_by_func(entry2, gui_zmatrix_line_changed, dialog);
g_signal_handlers_unblock_by_func(entry3, gui_zmatrix_line_changed, dialog);
}

/*********************************************/
/* compute the zmatrix for selected molecule */
/*********************************************/
void gui_zmatrix_compute(GtkWidget *w, gpointer dialog)
{
zmat_build();
gui_zmatrix_text_update(NULL, dialog);
gui_zmat_line_update(NULL, dialog);
}

/*************************************/
/* zmatrix coordinate editing widget */
/*************************************/
void gui_zmatrix_widget(GtkWidget *box, gpointer dialog)
{
GList *keys=NULL;
GtkWidget *swin, *view, *hbox, *label, *spin, *entry1, *entry2, *entry3;
GtkTextBuffer *buffer;
struct zmat_pak *zmat;
struct model_pak *model = dialog_model(dialog);

g_assert(model != NULL);
zmat = model->zmatrix;
g_assert(zmat != NULL);

/* scrolled text window for zmatrix lines */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
view = gtk_text_view_new();
gtk_text_view_set_editable(GTK_TEXT_VIEW(view), FALSE);
gtk_container_add(GTK_CONTAINER(swin), view);
gtk_box_pack_start(GTK_BOX(box), swin, TRUE, TRUE, 0);
buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
gtk_text_buffer_create_tag(buffer, "bg_changed", "background", "tan", NULL); 
dialog_child_set(dialog, "view1", view);
gui_zmatrix_text_update(NULL, dialog);


/* normal zmatrix line editing entries */
hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);

spin = gtk_spin_button_new_with_range(1, 999, 1);
gtk_box_pack_start(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
dialog_child_set(dialog, "spin1", spin);

entry1 = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry1), "blah");
gtk_box_pack_start(GTK_BOX(hbox), entry1, FALSE, FALSE, 0);
dialog_child_set(dialog, "entry1", entry1);

entry2 = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry2), "blah");
gtk_box_pack_start(GTK_BOX(hbox), entry2, FALSE, FALSE, 0);
dialog_child_set(dialog, "entry2", entry2);

entry3 = gtk_entry_new();
gtk_entry_set_text(GTK_ENTRY(entry3), "blah");
gtk_box_pack_start(GTK_BOX(hbox), entry3, FALSE, FALSE, 0);
dialog_child_set(dialog, "entry3", entry3);

gui_zmat_line_update(NULL, dialog);

g_signal_connect(GTK_OBJECT(spin), "value-changed", GTK_SIGNAL_FUNC(gui_zmat_line_update), dialog);

g_signal_connect(GTK_OBJECT(entry1), "changed",
                 GTK_SIGNAL_FUNC(gui_zmatrix_line_changed), dialog);
g_signal_connect(GTK_OBJECT(entry2), "changed",
                 GTK_SIGNAL_FUNC(gui_zmatrix_line_changed), dialog);
g_signal_connect(GTK_OBJECT(entry3), "changed",
                 GTK_SIGNAL_FUNC(gui_zmatrix_line_changed), dialog);


/* zmatrix variable editing section */
g_hash_table_foreach(zmat->vars, &gui_zmatrix_vars, &keys);
if (keys)
  {
  hbox = gtk_hbox_new(FALSE,0);
  gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
  label = gtk_label_new("Variables");
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

  hbox = gtk_hbox_new(FALSE,0);
  gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);

  zmat_combo = gtk_combo_new();
  gtk_combo_set_popdown_strings(GTK_COMBO(zmat_combo), keys);
  gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(zmat_combo)->entry), FALSE);
  gtk_box_pack_start(GTK_BOX(hbox), zmat_combo, FALSE, FALSE, 0);

  g_signal_connect(GTK_OBJECT(GTK_COMBO(zmat_combo)->entry), "changed",
                   (gpointer) gui_zmatrix_value_change, model->zmatrix);

  zmat_entry = gtk_entry_new();
  gtk_box_pack_start(GTK_BOX(hbox), zmat_entry, FALSE, FALSE, 0);

  g_signal_connect(GTK_OBJECT(zmat_entry), "changed",
                   (gpointer) gui_zmatrix_value_update, model->zmatrix);

  gui_zmatrix_value_change(NULL, model->zmatrix);
  }

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, PANEL_SPACING);
gui_button_x("Recompute geometry ", gui_zmatrix_model_update, dialog, hbox);

hbox = gtk_hbox_new(FALSE,0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
gui_button_x("Build zmatrix from selection ", gui_zmatrix_compute, dialog, hbox);
}

/*************************************/
/* zmatrix coordinate editing dialog */
/*************************************/
/* TODO - incorporate in terry's dialog or general coord editing */
void gui_zmat_dialog(void)
{
gpointer dialog;
GtkWidget *window, *vbox, *hbox, *label;
struct model_pak *model;
struct zmat_pak *zmat;

model = sysenv.active_model;
if (!model)
  return;

zmat = model->zmatrix;
if (!zmat)
  {
/* TODO - create? */
  printf("ERROR: empty zmatrix\n");
  return;
  }

/* request an energetics dialog */
dialog = dialog_request(ZMATRIX, "ZMATRIX editor", NULL, NULL, model);
if (!dialog)
  return;
window = dialog_window(dialog);

gtk_window_set_default_size(GTK_WINDOW(window), -1, 600);


/* units information */
vbox = gui_frame_vbox(NULL, FALSE, FALSE, GTK_DIALOG(window)->vbox);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

label = gtk_label_new("Distance units: ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
label = gtk_label_new(zmat_distance_units_get(model->zmatrix));
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

label = gtk_label_new("Angle units: ");
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
label = gtk_label_new(zmat_angle_units_get(model->zmatrix));
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);


/* main display */
vbox = gui_frame_vbox(NULL, TRUE, TRUE, GTK_DIALOG(window)->vbox);

gui_zmatrix_widget(vbox, dialog);

gtk_widget_show_all(window);
}
