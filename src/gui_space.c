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
#include <strings.h>
#include <ctype.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "space.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

/* main structure */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/************************/
/* draw periodic images */
/************************/
void gtk_refresh_images(GtkWidget *w, gpointer dummy)
{
struct model_pak *model;

model = sysenv.active_model;
g_assert(model != NULL);

/* atom colour/display updates */
space_make_images(CREATE, model);
coords_init(CENT_COORDS, model);
redraw_canvas(SINGLE);
}

/****************************/
/* display property updates */
/****************************/
void cb_make_asymmetric(GtkWidget *w, gpointer dummy)
{
GSList *list;
struct core_pak *core;
struct shel_pak *shel;
struct model_pak *model;

model = sysenv.active_model;

g_assert(model != NULL);

if (model->periodic)
  {
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    if (core->primary)
      continue;

    if (model->asym_on)
      core->status |= HIDDEN;
    else
      core->status &= ~HIDDEN;
    }
  for (list=model->shels ; list ; list=g_slist_next(list))
    {
    shel = list->data;
    if (shel->primary)
      continue;
    if (model->asym_on)
      shel->status |= HIDDEN;
    else
      shel->status &= ~HIDDEN;
    }
  }

/* update */
coords_compute(model);
redraw_canvas(SINGLE);
}

/***************************/
/* perioidic image globals */
/***************************/
gboolean image_update = TRUE;
GtkWidget *image_spin[6];

/****************************/
/* periodic image callbacks */
/****************************/
void cb_space_image_spinner(GtkWidget *w, gpointer data)
{
gint i;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return;

/* alter the spinner value */
i = GPOINTER_TO_INT(data);
g_assert(i >= 0);
g_assert(i < 6);
model->image_limit[i] = SPIN_FVAL(GTK_SPIN_BUTTON(w));

/* update only if this is a genuine call for perioidic image creation */
if (image_update)
  {
  space_make_images(CREATE, model);
  coords_init(CENT_COORDS, model);
  redraw_canvas(SINGLE);
  }
}

/********************************/
/* periodic image widget update */
/********************************/
void space_image_widget_redraw(void)
{
gint i, periodic;
struct model_pak *model;

/* setup */
model = sysenv.active_model;
if (!model)
  periodic = 0;
else
  periodic = model->periodic;

/* NB: don't allow updates as we're not making periodic images, */
/* we're setting the widget values to what they should be */
image_update = FALSE;
for (i=0 ; i<3 ; i++)
  {
  if (i < periodic)
    {
/* show and set available spinners */
    gtk_widget_set_sensitive(GTK_WIDGET(image_spin[2*i]), TRUE);
    gtk_widget_set_sensitive(GTK_WIDGET(image_spin[2*i+1]), TRUE);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(image_spin[2*i]), model->image_limit[2*i]);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(image_spin[2*i+1]), model->image_limit[2*i+1]);
    }
  else
    {
/* reset unavailable spinners to default values */
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(image_spin[2*i]), 0.0);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(image_spin[2*i+1]), 1.0);
/* hide unavailable spinners */
    gtk_widget_set_sensitive(GTK_WIDGET(image_spin[2*i]), FALSE);
    gtk_widget_set_sensitive(GTK_WIDGET(image_spin[2*i+1]), FALSE);
    }
  }
image_update = TRUE;
}

/**********************************/
/* periodic image default setting */
/**********************************/
void space_image_widget_reset(void)
{
struct model_pak *model;

/* update model (if exists) */
model = sysenv.active_model;
if (model)
  {
  space_make_images(INITIAL, model);
  coords_init(CENT_COORDS, model);
  }

/* update widget and canvas */
gui_refresh(GUI_MODEL_PROPERTIES);
gui_refresh(GUI_CANVAS);
}

/*********************************************/
/* periodic image creation for the main pane */
/*********************************************/
void space_image_widget_setup(GtkWidget *box)
{
gint i;
GtkWidget *table, *label;

/* use table for nice spacing */
table = gtk_table_new(3, 4, FALSE);
gtk_box_pack_start(GTK_BOX(box), table, FALSE, FALSE, 0);

/* image directions */
label = gtk_label_new ("-");
gtk_table_attach_defaults(GTK_TABLE(table),label,1,2,0,1);
label = gtk_label_new ("+");
gtk_table_attach_defaults(GTK_TABLE(table),label,2,3,0,1);

/* labels and spinners */
for (i=0 ; i<3 ; i++)
  {
/* negative direction */
  image_spin[2*i] = new_spinner(NULL, 0, 10, 1, 
                                cb_space_image_spinner,
                                GINT_TO_POINTER(2*i),
                                NULL);
/* TODO - set these to 1, when fractional periodic images are working */
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(image_spin[2*i]), 0);
  gtk_table_attach_defaults(GTK_TABLE(table), image_spin[2*i], 1, 2, i+1, i+2);

/* positive direction */
  image_spin[2*i+1] = new_spinner(NULL, 1, 10, 1,
                                  cb_space_image_spinner,
                                  GINT_TO_POINTER(2*i+1),
                                  NULL);
/* TODO - set these to 1, when fractional periodic images are working */
  gtk_spin_button_set_digits(GTK_SPIN_BUTTON(image_spin[2*i+1]), 0);
  gtk_table_attach_defaults(GTK_TABLE(table),
                            image_spin[2*i+1], 2, 3, i+1, i+2);
  }

/* set the spinner values */
gui_refresh(GUI_MODEL_PROPERTIES);
}

/***************************/
/* Space Group info dialog */
/***************************/
void periodicity_page(GtkWidget *box)
{
GString *info, *item;
GtkWidget *frame, *vbox;
struct model_pak *model;

/* checks */
model = sysenv.active_model;
if (!model)
  return;
if (!model->periodic)
  return;

item = g_string_new(NULL);
info = g_string_new(NULL);

#if GENPOS
/* FRAME */
frame = gtk_frame_new ("General positions");
gtk_box_pack_start(GTK_BOX(main_box),frame,TRUE,TRUE,0);
gtk_container_set_border_width (GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new (FALSE, 5);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
gtk_container_add (GTK_CONTAINER(frame),vbox);
/* table for nice spacing */
/* more intelligent dimension choice? */
cols = 4;
rows = model->sginfo.order / cols;
table = gtk_table_new(cols, rows, FALSE);
gtk_container_add(GTK_CONTAINER(vbox), table);

/* order of the group (ie mth matrix) */
m=0;
/* loop over rows */
while (m<model->sginfo.order)
  {
/* get row(i) and column(j) */
  i = m / cols;
  j = m % cols;
/* construct the three generator elements */
  mat = *(model->sginfo.matrix+m);
  g_string_sprintf(info," ");
  for (k=0 ; k<3 ; k++)
    {
/* offset elements */
    value = *(*(model->sginfo.offset+m)+k);
    need_sign = 0;
/* separate 2nd and 3rd elements with a comma */
    if (k)
      g_string_append(info, ", ");
/* clear the string for the new element to append */
    g_string_assign(item, "");
/* append the offset */
    if (fabs(2.0*value-1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "1/2");
      need_sign++;
      }
    if (fabs(2.0*value+1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "-1/2");
      need_sign++;
      }
    if (fabs(3.0*value-1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "1/3");
      need_sign++;
      }
    if (fabs(3.0*value+1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "-1/3");
      need_sign++;
      }
    if (fabs(4.0*value-1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "1/4");
      need_sign++;
      }
    if (fabs(4.0*value+1.0) < FRACTION_TOLERANCE)
      {
      g_string_sprintfa(item, "-1/4");
      need_sign++;
      }
/* append the appropriate matrix element(s) */
    for (l='x' ; l<='z' ; l++)
      {
      if (*mat > FRACTION_TOLERANCE)
        {
        if (need_sign)
          g_string_sprintfa(item, "+%c", l);
        else
          {
          g_string_sprintfa(item, "%c", l);
          need_sign++;
          }
        }
      if (*mat < -FRACTION_TOLERANCE)
        {
        g_string_sprintfa(item, "-%c",l);
        need_sign++;
        }
      mat++;
      }
    g_string_append(info, item->str);
    }
/* write the generator to the table */
  label = gtk_label_new(info->str);
  gtk_table_attach_defaults(GTK_TABLE(table),label,j,j+1,i,i+1);
  m++;
  }
#endif

/* FRAME */
frame = gtk_frame_new(NULL);
gtk_box_pack_start(GTK_BOX (box), frame, FALSE, FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
vbox = gtk_vbox_new(TRUE, PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
gtk_container_set_border_width(GTK_CONTAINER(GTK_BOX(vbox)), PANEL_SPACING);

gui_auto_check("Asymmetric cell", cb_make_asymmetric, NULL, &model->asym_on, vbox);

/* done, draw, exit */
gtk_widget_show_all(box);

g_string_free(info, TRUE);
g_string_free(item, TRUE);
}

