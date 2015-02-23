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
#include <math.h>

#include "gdis.h"
#include "dialog.h"
#include "interface.h"
#include "gui_image.h"

#include "go.xpm"
#include "pause.xpm"
#include "play.xpm"
#include "rewind.xpm"
#include "fastforward.xpm"
#include "stop.xpm"
#include "step_forward.xpm"
#include "step_backward.xpm"

extern struct sysenv_pak sysenv;

extern GtkWidget *window;

enum {AUTO_CHECK, AUTO_SPIN, AUTO_RANGE,
      AUTO_TEXT_ENTRY, AUTO_TEXT_LABEL,
      AUTO_INT_LABEL, AUTO_FLOAT_LABEL};

struct relation_pak 
{
/* true  - direct variable <-> widget relation */
/* false - relative offset (ie relation depends on active model) */
gboolean direct;

/* related object data */
struct model_pak *model;
gpointer variable;
gint offset;
GtkWidget *widget;
/* widget type (for updating) ie spin, check */
/* TODO - eliminate this using if GTK_IS_SPIN_BUTTON/CHECK_BUTTON macros */
gint type;
};

GSList *gui_relation_list=NULL;

#define DEBUG_RELATION 0

/**************************************/
/* auto update widget handling events */
/**************************************/
void gui_relation_destroy(GtkWidget *w, gpointer dummy)
{
GSList *list;
struct relation_pak *relation;

#if DEBUG_RELATION
printf("destroying widget (%p)", w);
#endif

/* find appropriate relation */
list = gui_relation_list;
while (list)
  {
  relation = list->data;
  list = g_slist_next(list);

/* free associated data */
  if (relation->widget == w)
    {
#if DEBUG_RELATION
printf(" : relation (%p)", relation);
#endif
    g_free(relation);
    gui_relation_list = g_slist_remove(gui_relation_list, relation);
    }
  }
#if DEBUG_RELATION
printf("\n");
#endif
}

/***************************************************/
/* set a value according to supplied widget status */
/***************************************************/
void gui_relation_set_value(GtkWidget *w, struct model_pak *model)
{
gpointer value;
GSList *list;
struct relation_pak *relation;

g_assert(w != NULL);

/* if no model supplied, use active model */
if (!model)
  model = sysenv.active_model;

/* find appropriate relation (direct -> widget, else model) */
list = gui_relation_list;
while (list)
  {
  relation = list->data;
  list = g_slist_next(list);

  if (relation->widget != w)
    continue;

  if (relation->direct)
    value = relation->variable;
  else
    {
    g_assert(model != NULL);
    value = (gpointer) model + relation->offset;
    }

/* update variable associated with the widget */
  switch (relation->type)
    {
    case AUTO_CHECK:
#if DEBUG_RELATION
printf("model %p : relation %p : setting variable to %d\n", 
       model, relation, gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)));
#endif
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
        *((gint *) value) = TRUE;
      else
        *((gint *) value) = FALSE;
      break;

    case AUTO_SPIN:
#if DEBUG_RELATION
printf("model %p : relation %p : setting variable to %f\n", 
       model, relation, SPIN_FVAL(GTK_SPIN_BUTTON(w)));
#endif
      *((gdouble *) value) = SPIN_FVAL(GTK_SPIN_BUTTON(w));
      break;

    case AUTO_RANGE:
      *((gint *) value) = gtk_range_get_value(GTK_RANGE(w));
      break;

    case AUTO_TEXT_ENTRY:
/* CURRENT */
      g_free(*((gchar **) value));
      *((gchar **) value) = g_strdup(gtk_entry_get_text(GTK_ENTRY(w)));
      break;

    }
  }
}

/***************************************************/
/* updates dependant widget with new variable data */
/***************************************************/
void gui_relation_update_widget(gpointer value)
{
gchar *text;
GSList *list;
struct relation_pak *relation;

/* loop over all existing relations */
for (list=gui_relation_list ; list ; list=g_slist_next(list))
  {
  relation = list->data;

  g_assert(relation != NULL);

  if (!relation->direct)
    continue;
  if (relation->variable != value)
    continue;

#if DEBUG_RELATION
printf("relation %p : type %d : updating widget\n", relation, relation->type);
#endif

/* synchronize widget and variable */
  switch (relation->type)
    {
    case AUTO_CHECK:
      if (*((gint *) value))
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(relation->widget),
                                     TRUE);
      else
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(relation->widget),
                                     FALSE);
      break;

    case AUTO_SPIN:
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(relation->widget),
                                *((gdouble *) value));
      break;

    case AUTO_RANGE:
      gtk_range_set_value(GTK_RANGE(relation->widget), *((gint *) value));
      break;

    case AUTO_TEXT_ENTRY:
/* CURRENT */
      text = *((gchar **) value);
      if (text)
        gtk_entry_set_text(GTK_ENTRY(relation->widget), text);
      break;

    case AUTO_TEXT_LABEL:
      gtk_label_set_text(GTK_LABEL(relation->widget), *((gchar **) value));
      break;

    case AUTO_INT_LABEL:
      text = g_strdup_printf("%d", *((gint *) value));
      gtk_label_set_text(GTK_LABEL(relation->widget), text);
      g_free(text);
      break;

    case AUTO_FLOAT_LABEL:
      text = g_strdup_printf("%.4f", *((gdouble *) value));
      gtk_label_set_text(GTK_LABEL(relation->widget), text);
      g_free(text);
      break;
    }
  }
}

/*****************************/
/* called for a model switch */
/*****************************/
/* ie sets all widget states according to related variable */
/* NB: pass NULL to get everything updated */
void gui_relation_update(gpointer data)
{
gchar *text;
GSList *list;
struct relation_pak *relation;
gpointer value;

/* loop over all existing relations */
for (list=gui_relation_list ; list ; list=g_slist_next(list))
  {
  relation = list->data;

  g_assert(relation != NULL);

/* update test */
  if (relation->direct)
    {
    if (data)
      if (data != relation->model)
        continue;
    value = relation->variable;
    }
  else
    {
    if (data)
      value = data + relation->offset;
    else
      continue;
    }

#if DEBUG_RELATION
printf("relation %p : type %d : model change\n", relation, relation->type);
#endif

/* synchronize widget and variable */
  switch (relation->type)
    {
    case AUTO_CHECK:
      if (*((gint *) value))
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(relation->widget),
                                     TRUE);
      else
        gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(relation->widget),
                                     FALSE);
      break;

    case AUTO_SPIN:
      gtk_spin_button_set_value(GTK_SPIN_BUTTON(relation->widget),
                                *((gdouble *) value));
      break;

    case AUTO_RANGE:
      gtk_range_set_value(GTK_RANGE(relation->widget), *((gint *) value));
      break;

    case AUTO_TEXT_ENTRY:
/* CURRENT */
      text = *((gchar **) value);
      if (text)
        gtk_entry_set_text(GTK_ENTRY(relation->widget), text);
      break;

    case AUTO_TEXT_LABEL:
      gtk_label_set_text(GTK_LABEL(relation->widget), *((gchar **) value));
      break;

    case AUTO_INT_LABEL:
      text = g_strdup_printf("%d", *((gint *) value));
      gtk_label_set_text(GTK_LABEL(relation->widget), text);
      g_free(text);
      break;

    case AUTO_FLOAT_LABEL:
      text = g_strdup_printf("%.4f", *((gdouble *) value));
      gtk_label_set_text(GTK_LABEL(relation->widget), text);
      g_free(text);
      break;
    }
  }
}

/**************************************************/
/* create a new variable to widget correspondance */
/**************************************************/
void
gui_relation_submit(GtkWidget        *widget,
                      gint              type,
                      gpointer          value,
                      gboolean          direct,
                      struct model_pak *model)
{
struct relation_pak *relation;

/* alloc and init relation */
relation = g_malloc(sizeof(struct relation_pak));
relation->direct = direct;
relation->type = type;
if (model)
  relation->offset = value - (gpointer) model;
else
  relation->offset = 0;
relation->variable = value;
relation->model = model;
relation->widget = widget;

/* this can happen if gui_auto_check() is used on a variable that is not in a model */
/* eg in sysenv ... in this case use gui_direct_check() */
/* FIXME - a better way to fail gracefully? better still - rewrite old code */
if (!direct)
  g_assert(relation->offset > 0);

#if DEBUG_RELATION
printf("submit: ");
printf("[%p type %d]", relation, type);
if (model)
  {
  if (direct)
    printf("[model %d][value %p]", model->number, value);
  else
    printf("[model %d][value %p][offset %d]", model->number, value, relation->offset);
  }
printf("\n");
#endif

switch (relation->type)
  {
  case AUTO_CHECK:
    if (*((gint *) value))
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), TRUE);
    else
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget), FALSE);
    break;

  case AUTO_SPIN:
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(widget), *((gdouble *) value));
    break;

  case AUTO_RANGE:
    gtk_range_set_value(GTK_RANGE(widget), *((gint *) value));
    break;
  }

gui_relation_list = g_slist_prepend(gui_relation_list, relation);
}

/*****************************************/
/* convenience routine for check buttons */
/*****************************************/
GtkWidget *new_check_button(gchar *label, gpointer callback, gpointer arg,
                                            gint state, GtkWidget *box)
{
GtkWidget *button;

/* create, pack & show the button */
button = gtk_check_button_new_with_label(label);
gtk_box_pack_start(GTK_BOX(box), button, TRUE, TRUE, 0);
gtk_widget_show(button);

/* set the state (NB: before the callback is attached) */
if (state)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

/* attach the callback */
if (callback)
  g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(callback), arg);
return(button);
}

/*********************************/
/* relative updated check button */
/*********************************/
/* NB: assumes gui_mode_switchl() will call gui_relation_update() */
GtkWidget *gui_auto_check(gchar *label, gpointer cb, gpointer arg,
                            gint *state, GtkWidget *box)
{
GtkWidget *button;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

/* create the button */
button = gtk_check_button_new_with_label(label);
gui_relation_submit(button, AUTO_CHECK, state, FALSE, model);
if (box)
  gtk_box_pack_start(GTK_BOX(box), button, TRUE, TRUE, 0);

/* callback to set the variable to match the widget */
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(gui_relation_set_value), NULL);

/* callback to do (user defined) update tasks */
if (cb)
  g_signal_connect_after(GTK_OBJECT(button), "clicked", cb, arg);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(button), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(button);
}

/*************************************************/
/* activate/deactivate based on a checkbox state */
/*************************************************/
/* TODO - shortcut for setup up one of these widgets */
void gui_checkbox_refresh(GtkWidget *w, GtkWidget *box)
{
if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w)))
  gtk_widget_set_sensitive(box, TRUE);
else
  gtk_widget_set_sensitive(box, FALSE);
}

/*********************************/
/* directly updated check button */
/*********************************/
/* NB: assumes gui_mode_switch() will call gui_relation_update() */
GtkWidget *gui_direct_check(gchar *label, gint *state,
                              gpointer cb, gpointer arg,
                              GtkWidget *box)
{
GtkWidget *button;
struct model_pak *model;

model = sysenv.active_model;

/* create the button */
button = gtk_check_button_new_with_label(label);
gui_relation_submit(button, AUTO_CHECK, state, TRUE, model);
if (box)
  gtk_box_pack_start(GTK_BOX(box), button, TRUE, TRUE, 0);


/* callback to set the variable to match the widget */
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(gui_relation_set_value), model);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(button), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

/* callback to do (user defined) update tasks */
if (cb)
  g_signal_connect_after(GTK_OBJECT(button), "clicked", cb, arg);

return(button);
}

/****************************/
/* a simple labelled button */
/****************************/
GtkWidget *gui_button(gchar *txt, gpointer cb, gpointer arg, GtkWidget *w, gint mask)
{
gint fill1, fill2;
GtkWidget *button;

/* create the button */
button = gtk_button_new_with_label(txt);

/* setup the packing requirements */
if (w)
  {
  fill1 = fill2 = FALSE;
  if (mask & 1)
    fill1 = TRUE;
  if (mask & 2)
    fill2 = TRUE;
  if (mask & 4)
    gtk_box_pack_end(GTK_BOX(w), button, fill1, fill2, 0);
  else
    gtk_box_pack_start(GTK_BOX(w), button, fill1, fill2, 0);
  }

/* attach the callback */
if (cb)
  g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(cb), arg);

return(button);
}

/************************************/
/* text label plus an X push button */
/************************************/
GtkWidget *gui_button_x(gchar *text,
                          gpointer cb, gpointer arg,
                          GtkWidget *box)
{
GtkWidget *hbox, *button, *label, *image;
GdkPixmap *pixmap;
GdkBitmap *mask;
GtkStyle *style;

/* create button */
button = gtk_button_new();
hbox = gtk_hbox_new(FALSE, 0);
gtk_container_add(GTK_CONTAINER(button), hbox);

/* create image */
style = gtk_widget_get_style(window);
pixmap = gdk_pixmap_create_from_xpm_d
          (window->window, &mask, &style->white, go_xpm);
image = gtk_image_new_from_pixmap(pixmap, mask);

gtk_box_pack_start(GTK_BOX(hbox), image, FALSE, FALSE, 0);

/* create label */
if (box)
  {
/* packing sub-widget */
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);
  if (text)
    {
    label = gtk_label_new(text);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    }
  gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
  }

if (cb)
  g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(cb), arg);

return(button);
}

/****************************************/
/* text label plus labelled push button */
/****************************************/
void gui_button_label(gchar *t1, gchar *t2, gpointer cb, gpointer arg, GtkWidget *box)
{
GtkWidget *hbox, *button, *label;

/* packing sub-widget */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);

/* create label */
label = gtk_label_new(t1);
gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);

/* create button */
button = gtk_button_new_with_label(t2);
g_signal_connect_swapped(GTK_OBJECT(button), "clicked", 
                         GTK_SIGNAL_FUNC(cb), arg);
gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
}

/*****************************************/
/* convenience routine for radio buttons */
/*****************************************/
gint active=0, count=0, fill1=FALSE, fill2=FALSE;
GSList *group=NULL;
GtkWidget *box=NULL;

void new_radio_group(gint i, GtkWidget *pack_box, gint mask)
{
group = NULL;
active = i;
count = 0;
box = pack_box;

/* setup the packing requirements */
fill1 = fill2 = FALSE;
if (mask & 1)
  fill1 = TRUE;
if (mask & 2)
  fill2 = TRUE;
}

GtkWidget *add_radio_button(gchar *label, gpointer call, gpointer arg)
{
GtkWidget *button;

/*
printf("Adding button: %d (%d) to (%p)\n", count, active, group);
*/

/* better error checking - even call new_radio_group() ??? */
g_return_val_if_fail(box != NULL, NULL);

/* make a new button */
button = gtk_radio_button_new_with_label(group, label);

/* get the group (for the next button) */
/* NB: it is vital that this is ALWAYS done */
group = gtk_radio_button_group(GTK_RADIO_BUTTON(button));

/* attach the callback */
g_signal_connect_swapped(GTK_OBJECT(button), "pressed",
                         GTK_SIGNAL_FUNC(call), (gpointer) arg);

/* pack the button */
gtk_box_pack_start(GTK_BOX(box), button, fill1, fill2, 0);

/* set active? */
if (active == count++)
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);

return(button);
}

/************************************/
/* create a colour selection dialog */
/************************************/
GtkWidget *new_csd(gchar *title, gpointer callback)
{
GtkWidget *csd;

/* create colour selection dialog */
csd = gtk_color_selection_dialog_new(title);

/* setup callbacks */
g_signal_connect_swapped(G_OBJECT(GTK_COLOR_SELECTION_DIALOG(csd)->ok_button),
                        "clicked", 
                         GTK_SIGNAL_FUNC(callback),
                        (gpointer) GTK_COLOR_SELECTION_DIALOG(csd)->colorsel);

g_signal_connect_swapped(G_OBJECT(GTK_COLOR_SELECTION_DIALOG(csd)->cancel_button),
                        "clicked",
                         GTK_SIGNAL_FUNC(gtk_widget_destroy),
                        (gpointer) csd);

/* done */
gtk_widget_show(csd);
return(csd);
}

/******************************/
/* create a label + a spinner */
/******************************/
GtkWidget *new_spinner(gchar *text, gdouble min, gdouble max, gdouble step,
                       gpointer callback, gpointer arg, GtkWidget *box)
{
GtkWidget *hbox, *label, *spin;

spin = gtk_spin_button_new_with_range(min, max, step);

/* optional */
if (box)
  {
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);
  if (text)
    {
    label = gtk_label_new(text);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    }
  gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
  }

if (callback)
  g_signal_connect_after(GTK_OBJECT(spin), "value-changed",
                         GTK_SIGNAL_FUNC(callback), arg);

return(spin);
}

/************************************************************/
/* as above, but automatically attaches spinner to callback */
/************************************************************/
GtkWidget *gui_new_spin(gchar *text, gdouble x0, gdouble x1, gdouble dx,
                          gpointer callback, GtkWidget *box)
{
GtkWidget *hbox, *label, *spin;

spin = gtk_spin_button_new_with_range(x0, x1, dx);
gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin), x0);

/* optional */
if (box)
  {
/* create the text/spinner layout */
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);
  if (text)
    {
    label = gtk_label_new(text);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
    }
  gtk_box_pack_end(GTK_BOX(hbox), spin, FALSE, FALSE, 0);
  }
if (callback)
  g_signal_connect_after(GTK_OBJECT(spin), "value-changed",
                         GTK_SIGNAL_FUNC(callback), spin);
return(spin);
}

/*********************************/
/* automatically updated spinner */
/*********************************/
/* current model relative correlation between value and spinner */
/* NB: assumes gui_mode_switchl() will call gui_relation_update() */
GtkWidget *gui_auto_spin(gchar *text, gdouble *value,
                           gdouble min, gdouble max, gdouble step,
                           gpointer callback, gpointer arg, GtkWidget *box)
{
GtkWidget *spin;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

/* create the text/spinner combo */
spin = new_spinner(text, min, max, step, callback, arg, box);

/* set up a relationship */
gui_relation_submit(spin, AUTO_SPIN, value, FALSE, model);

/* callback to set the variable to match the widget */
g_signal_connect(GTK_OBJECT(spin), "value-changed",
                (gpointer) gui_relation_set_value, NULL);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(spin), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(spin);
}

/*********************************/
/* automatically updated spinner */
/*********************************/
/* direct correlation between value and spinner */
GtkWidget *gui_direct_spin(gchar *text, gdouble *value,
                             gdouble min, gdouble max, gdouble step,
                             gpointer callback, gpointer arg, GtkWidget *box)
{
gint size=0;
gdouble digits;
GtkWidget *spin;
struct model_pak *model;

model = sysenv.active_model;

/* create the text/spinner combo */
spin = new_spinner(text, min, max, step, NULL, NULL, box);

/* HACK - cope with GTK underestimating the size needed to display spin values */
/* TODO - examine sig fig of *value to get dp */
digits = log10(step);
if (digits < 0.0)
  {
  digits -= 0.9;
  size = 3 + fabs(digits);
  }
if (size < 5)
  size = 5;

gtk_widget_set_size_request(spin, sysenv.gtk_fontsize*size, -1);

/*
printf("%f : %f -> %f (%d)\n", max, step, digits, size);
gtk_widget_set_size_request(spin, sysenv.gtk_fontsize*5, -1);
*/

/* set up a relationship */
gui_relation_submit(spin, AUTO_SPIN, value, TRUE, model);

/* callback to set the variable to match the widget */
g_signal_connect(GTK_OBJECT(spin), "value-changed",
                 GTK_SIGNAL_FUNC(gui_relation_set_value), model);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(spin), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

/* connect after, so all updates done before user callback is invoked */
if (callback)
  g_signal_connect_after(GTK_OBJECT(spin), "value-changed",
                         GTK_SIGNAL_FUNC(callback), arg);

return(spin);
}

/**************************/
/* auto update slider bar */
/**************************/
GtkWidget *gui_direct_hscale(gdouble min, gdouble max, gdouble step,
                               gpointer value, gpointer func, gpointer arg,
                                                            GtkWidget *box)
{
GtkWidget *hscale;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

hscale = gtk_hscale_new_with_range(min, max, step);
gtk_range_set_update_policy(GTK_RANGE(hscale), GTK_UPDATE_CONTINUOUS);

if (box)
  gtk_box_pack_start(GTK_BOX(box), hscale, TRUE, TRUE, 0);

gui_relation_submit(hscale, AUTO_RANGE, value, TRUE, model);

/* callback to set the variable to match the widget */
g_signal_connect(GTK_OBJECT(hscale), "value_changed",
                 GTK_SIGNAL_FUNC(gui_relation_set_value), model);

/* user callback */
if (func)
  g_signal_connect_after(GTK_OBJECT(hscale), "value_changed",
                         GTK_SIGNAL_FUNC(func), arg);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(hscale), "destroy",
                GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(hscale);
}

/**************************************************/
/* convenience function for labelled text entries */
/**************************************************/
GtkWidget *gui_text_entry(gchar *text, gchar **value,
                          gint edit, gint pack, GtkWidget *box)
{
GtkWidget *hbox, *label, *entry;

g_assert(box != NULL);

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
if (text)
  {
  label = gtk_label_new(text);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  }
entry = gtk_entry_new();

if (value)
  {
  if (*value)
    gtk_entry_set_text(GTK_ENTRY(entry), *value);
  }

gtk_entry_set_editable(GTK_ENTRY(entry), edit);

if (value)
  gui_relation_submit(entry, AUTO_TEXT_ENTRY, value, TRUE, NULL);

gtk_box_pack_end(GTK_BOX(hbox), entry, pack, TRUE, 0);

if (value)
  {
  if (edit)
    g_signal_connect(GTK_OBJECT(entry), "changed",
                     GTK_SIGNAL_FUNC(gui_relation_set_value), NULL);

  g_signal_connect(GTK_OBJECT(entry), "destroy",
                   GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);
  }

return(entry);
}

/**************************/
/* auto update text label */
/**************************/
GtkWidget *gui_auto_text_label(gchar **text)
{
GtkWidget *label;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);
g_assert(text != NULL);

if (*text)
  label = gtk_label_new(*text);
else
  label = gtk_label_new("null");

gui_relation_submit(label, AUTO_TEXT_LABEL, text, FALSE, model);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(label), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(label);
}

/**************************/
/* auto update text label */
/**************************/
GtkWidget *gui_auto_int_label(gint *value)
{
gchar *text;
GtkWidget *label;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

text = g_strdup_printf("%d", *value);
label = gtk_label_new(text);
g_free(text);

gui_relation_submit(label, AUTO_INT_LABEL, value, FALSE, model);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(label), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(label);
}

/**************************/
/* auto update text label */
/**************************/
GtkWidget *gui_auto_float_label(gdouble *value)
{
gchar *text;
GtkWidget *label;
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

text = g_strdup_printf("%.4f", *value);
label = gtk_label_new(text);
g_free(text);

gui_relation_submit(label, AUTO_FLOAT_LABEL, value, FALSE, model);

/* callback to remove the variable <-> widget relation */
g_signal_connect(GTK_OBJECT(label), "destroy",
                 GTK_SIGNAL_FUNC(gui_relation_destroy), NULL);

return(label);
}

/*****************************/
/* create a stock GTK button */
/*****************************/
GtkWidget * 
gui_stock_button(const gchar *id, gpointer cb, gpointer arg, GtkWidget *box)
{
GtkWidget *button;

button = gtk_button_new_from_stock(id);

if (box)
  gtk_box_pack_start(GTK_BOX(box), button, FALSE, FALSE, 0);

g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(cb), arg);

return(button);
}

/***************************************************/
/* create a button with a stock icon & custom text */
/***************************************************/
GtkWidget * 
gui_icon_button(const gchar *id, const gchar *text,
                  gpointer cb, gpointer arg,
                  GtkWidget *box)
{
gint stock = TRUE;
gpointer data;
GtkWidget *hbox, *button, *label, *image=NULL;
GdkBitmap *mask;
GdkPixmap *pixmap;
GtkStyle *style;

/* button */
button = gtk_button_new();
hbox = gtk_hbox_new(FALSE, 4);
gtk_container_add(GTK_CONTAINER(button), hbox);


/* CURRENT */
if (g_ascii_strncasecmp(id, "IMAGE_", 6) == 0)
  {
  data = image_table_lookup(id);
  if (data)
    {
    image = gtk_image_new_from_pixbuf(data);
    stock = FALSE;
    }
  else
    printf("[%s] : image not found.\n", id);
  }


/* GDIS icons */
if (g_ascii_strncasecmp(id, "GDIS_PAUSE", 10) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        pause_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_PLAY", 9) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        play_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_REWIND", 11) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        rewind_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_FASTFORWARD", 16) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        fastforward_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_STOP", 9) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        stop_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_STEP_FORWARD", 17) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        step_forward_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }
if (g_ascii_strncasecmp(id, "GDIS_STEP_BACKWARD", 18) == 0)
  {
  style = gtk_widget_get_style(window);

  pixmap = gdk_pixmap_create_from_xpm_d(window->window, &mask, &style->white,
                                        step_backward_xpm);
  image = gtk_image_new_from_pixmap(pixmap, mask);

  stock = FALSE;
  }

/* standard GTK */
if (stock)
  image = gtk_image_new_from_stock(id, GTK_ICON_SIZE_BUTTON);

/* label dependent packing  */
if (text)
  {
  gtk_box_pack_start(GTK_BOX(hbox), image, FALSE, FALSE, 0);
  label = gtk_label_new(text);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  }
else
  {
  gtk_box_pack_start(GTK_BOX(hbox), image, TRUE, FALSE, 0);
  }

/* packing & callback */
if (box)
  gtk_box_pack_start(GTK_BOX(box), button, TRUE, TRUE, 0);

if (cb)
  g_signal_connect(GTK_OBJECT(button), "clicked", GTK_SIGNAL_FUNC(cb), arg);

return(button);
}

/*******************************************/
/* update function for scolled text window */
/*******************************************/
void gui_text_window_update(GtkTextBuffer *buffer, gchar **text)
{
GtkTextIter start, end;

gtk_text_buffer_get_bounds(buffer, &start, &end);
g_free(*text);
*text = gtk_text_buffer_get_text(buffer, &start, &end, FALSE);
}

/*********************************/
/* short cut for displaying text */
/*********************************/
GtkWidget *gui_text_window(gchar **text, gint editable)
{
gint n;
GtkWidget *swin, *view;
GtkTextBuffer *buffer;

g_assert(text != NULL);

swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
view = gtk_text_view_new();
gtk_text_view_set_editable(GTK_TEXT_VIEW(view), editable);
gtk_container_add(GTK_CONTAINER(swin), view);
buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(view));
/* tags??? */
/*
gtk_text_buffer_create_tag(buffer, "fg_blue", "foreground", "blue", NULL);  
gtk_text_buffer_create_tag(buffer, "fg_red", "foreground", "red", NULL);
*/
if (*text)
  {
  n = strlen(*text);
  gtk_text_buffer_insert_at_cursor(buffer, *text, n);
  }

if (editable)
  g_signal_connect(G_OBJECT(buffer), "changed",
                   GTK_SIGNAL_FUNC(gui_text_window_update), text);

return(swin);
}

/**********************************/
/* vbox in a frame layout shorcut */
/**********************************/
GtkWidget *gui_frame_vbox(const gchar *label, gint p1, gint p2, GtkWidget *box)
{
GtkWidget *frame, *vbox;

g_assert(box != NULL);

frame = gtk_frame_new(label);
gtk_box_pack_start(GTK_BOX(box), frame, p1, p2, 0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), 0.5*PANEL_SPACING);
vbox = gtk_vbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(vbox), 0.5*PANEL_SPACING);

return(vbox);
}

/**********************************/
/* hbox in a frame layout shorcut */
/**********************************/
GtkWidget *gui_frame_hbox(const gchar *label, gint p1, gint p2, GtkWidget *box)
{
GtkWidget *frame, *hbox;

g_assert(box != NULL);

frame = gtk_frame_new(label);
gtk_box_pack_start(GTK_BOX(box), frame, p1, p2, 0); 
gtk_container_set_border_width(GTK_CONTAINER(frame), 0.5*PANEL_SPACING);
hbox = gtk_hbox_new(FALSE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), hbox);
gtk_container_set_border_width(GTK_CONTAINER(hbox), 0.5*PANEL_SPACING);

return(hbox);
}

/***************************************************************/
/* horizontal packing function for label plus some data widget */
/***************************************************************/
void gui_hbox_pack(GtkWidget *box, gchar *text, GtkWidget *w, gint pack)
{
gint pack1, pack2;
GtkWidget *label, *hbox;

g_assert(box != NULL);

/* TODO - the pack argument is intended to allow the possibility of different packing styles */
/* TODO - switch (pack) {} */
pack1 = TRUE;
pack2 = TRUE;

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, pack1, pack2, 0);

if (text)
  {
  label = gtk_label_new(text);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, TRUE, 0);
  gtk_box_pack_end(GTK_BOX(hbox), w, FALSE, FALSE, 0);
  }
else
  {
  if (w)
    gtk_box_pack_end(GTK_BOX(hbox), w, TRUE, TRUE, 0);
  }
}

/**************************************/
/* pulldown menu convenience function */
/**************************************/
/* TODO - when ready - will replace the old gui_pulldown_new() function */
gpointer gui_pd_new(GSList *list, gint active, gpointer callback, gpointer argument)
{
GSList *item;
GtkTreeModel *treemodel;
GtkTreePath *treepath;
GtkTreeIter iter;
GtkWidget *w = NULL;

#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 4
/* build the widget and associated drop down list */
w = gtk_combo_box_new_text();
for (item=list ; item ; item=g_slist_next(item))
  gtk_combo_box_append_text(GTK_COMBO_BOX(w), item->data);

/* set the currently active item */
treemodel = gtk_combo_box_get_model(GTK_COMBO_BOX(w));
treepath = gtk_tree_path_new_from_indices(active, -1);
gtk_tree_model_get_iter(treemodel, &iter, treepath);
gtk_combo_box_set_active_iter(GTK_COMBO_BOX(w), &iter);

/* attach callback (if any) */
if (callback)
  g_signal_connect(GTK_OBJECT(w), "changed", GTK_SIGNAL_FUNC(callback), argument);
#endif

return(w);
}

/*******************************************************/
/* retrieve the current active pulldown menu item text */
/*******************************************************/
gchar *gui_pd_text(gpointer pd)
{
#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 6
return(gtk_combo_box_get_active_text(GTK_COMBO_BOX(pd)));
#else
/* you're screwed */
return(NULL);
#endif
}

/**************************************/
/* pulldown menu convenience function */
/**************************************/
gpointer gui_pulldown_new(const gchar *text, GList *list, gint edit, GtkWidget *box)
{
GtkWidget *hbox, *label, *combo;

combo = gtk_combo_new();
gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);
gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), edit);

if (box)
  {
/* create the text/spinner layout */
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), hbox, TRUE, TRUE, 0);
  if (text)
    {
    label = gtk_label_new(text);
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, TRUE, 0);
    gtk_box_pack_end(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    }
  else
    gtk_box_pack_end(GTK_BOX(hbox), combo, TRUE, TRUE, 0);
  }
return(GTK_COMBO(combo)->entry);
}

/********************************/
/* pulldown menu data retrieval */
/********************************/
const gchar *gui_pulldown_text(gpointer w)
{
if (GTK_IS_ENTRY(w))
  return(gtk_entry_get_text(GTK_ENTRY(w)));
return(NULL);
}

/*******************************/
/* labelled colour editing box */
/*******************************/
void gui_colour_box(const gchar *text, gdouble *rgb, GtkWidget *box)
{
GtkWidget *hbox, *label, *button;
GdkColor colour;

hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);
if (text)
  {
  label = gtk_label_new(text);
  gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 0);
  }

button = gtk_button_new_with_label("    ");
gtk_box_pack_end(GTK_BOX(hbox), button, FALSE, FALSE, 0);
g_signal_connect(GTK_OBJECT(button), "clicked",
                 GTK_SIGNAL_FUNC(dialog_colour_new), rgb);

colour.red   = rgb[0]*65535.0;
colour.green = rgb[1]*65535.0;
colour.blue  = rgb[2]*65535.0;
gtk_widget_modify_bg(button, GTK_STATE_NORMAL, &colour);
}
