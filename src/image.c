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
#include <string.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "matrix.h"
#include "opengl.h"
#include "dialog.h"
#include "interface.h"
#include "gui_image.h"


#include "folder.xpm"
#include "disk.xpm"
#include "arrow.xpm"
#include "axes.xpm"
#include "tools.xpm"
#include "palette.xpm"
#include "cross.xpm"
#include "geom.xpm"
#include "cell.xpm"
#include "camera.xpm"
#include "element.xpm"
#include "tb_animate.xpm"
#include "tb_diffraction.xpm"
#include "tb_isosurface.xpm"
#include "tb_surface.xpm"
#include "canvas_single.xpm"
#include "canvas_create.xpm"
#include "canvas_delete.xpm"

/* NEW: include a button to create a new model */
#include "plus.xpm"
 


extern struct sysenv_pak sysenv;

extern GtkWidget *window;

/****************************************/
/* write an image of the current canvas */
/****************************************/
void image_write(gpointer w)
{
GdkPixbuf *pixbuf;
GError *error=NULL;

pixbuf = gdk_pixbuf_get_from_drawable(NULL, w, NULL,
                                      0, 0, 0, 0,
                                      sysenv.width, sysenv.height);

gdk_pixbuf_save(pixbuf, sysenv.snapshot_filename, "jpeg",
                &error, "quality", "90", NULL);

g_free(sysenv.snapshot_filename);
sysenv.snapshot = FALSE;
}

/****************************************/
/* callback to schedule a canvas export */
/****************************************/
void image_export(gchar *name)
{
g_assert(name != NULL);

dialog_destroy_type(FILE_SELECT);

/* FIXME - all GTK stuff (ie closing the file dialog/raising windows */
/* to the foreground etc.) will be done only AFTER this routine ends */
/* ie returns to the control of gtk_main_loop - so there is no way */
/* of preventing the dialog from getting in the way */
/* probably have to set a flag to do it after the next redraw_canvas() */
sysenv.snapshot = TRUE;
sysenv.snapshot_filename = g_build_filename(sysenv.cwd, name, NULL);
}

/*************************************************/
/* add a filename to active model's picture list */
/*************************************************/
void image_import(const gchar *name)
{
gchar *picture;
struct model_pak *model;

/* checks */
/* TODO - check file validity */
g_assert(name != NULL);

dialog_destroy_type(FILE_SELECT);

model = sysenv.active_model;
if (!model)
  return;

picture = g_strdup(name);

model->picture_list = g_slist_append(model->picture_list, picture);
model->picture_active = picture;
model->graph_active = NULL;

/* updates */
tree_model_add(model);
redraw_canvas(SINGLE);
}

/***********************************/
/* get a filename for image export */
/***********************************/
void image_export_dialog(void)
{
if (sysenv.active_model)
  file_dialog("Export image", NULL, FILE_SAVE, 
             (gpointer) image_export, PICTURE);
}

/*******************************************/
/* get a file dialog for picture selection */
/*******************************************/
void image_import_dialog(void)
{
if (sysenv.active_model)
  file_dialog("Import picture", NULL, FILE_LOAD, 
             (gpointer) image_import, PICTURE);
}

/**************************************/
/* retrieve a standard image's pixbuf */
/**************************************/
gpointer image_table_lookup(const gchar *name)
{
return(g_hash_table_lookup(sysenv.image_table, name));
}

/******************************************/
/* setup the internal pixbuf lookup table */
/******************************************/
void image_table_init(void)
{
gint i;
gchar *name;
GdkPixbuf *pixbuf;

sysenv.image_table = g_hash_table_new(g_str_hash, g_str_equal);

for (i=0 ; i<IMAGE_LAST ; i++)
  {
/* init the name and pixbuf data */
  switch (i)
    {
    case IMAGE_ANIMATE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) tb_animate_xpm);
      name = "image_animate";
      break;
    case IMAGE_ARROW:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) arrow_xpm);
      name = "image_arrow";
      break;
    case IMAGE_AXES:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) axes_xpm);
      name = "image_axes";
      break;
    case IMAGE_CANVAS_SINGLE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) canvas_single_xpm);
      name = "image_canvas_single";
      break;
    case IMAGE_CANVAS_CREATE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) canvas_create_xpm);
      name = "image_canvas_create";
      break;
    case IMAGE_CANVAS_DELETE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) canvas_delete_xpm);
      name = "image_canvas_delete";
      break;
    case IMAGE_CAMERA:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) camera_xpm);
      name = "image_camera";
      break;
    case IMAGE_COMPASS:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) geom_xpm);
      name = "image_compass";
      break;
    case IMAGE_CROSS:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) cross_xpm);
      name = "image_cross";
      break;
    case IMAGE_DIFFRACTION:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) tb_diffraction_xpm);
      name = "image_diffraction";
      break;
    case IMAGE_DISK:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) disk_xpm);
      name = "image_disk";
      break;
    case IMAGE_ELEMENT:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) element_xpm);
      name = "image_element";
      break;
    case IMAGE_FOLDER:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) folder_xpm);
      name = "image_folder";
      break;
    case IMAGE_ISOSURFACE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) tb_isosurface_xpm);
      name = "image_isosurface";
      break;
    case IMAGE_MEASURE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) geom_xpm);
      name = "image_measure";
      break;
    case IMAGE_PERIODIC:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) cell_xpm);
      name = "image_periodic";
      break;
    case IMAGE_PALETTE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) palette_xpm);
      name = "image_palette";
      break;
    case IMAGE_PLUS:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) plus_xpm);
      name = "image_plus";
      break;    
    case IMAGE_SURFACE:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) tb_surface_xpm);
      name = "image_surface";
      break;
    case IMAGE_TOOLS:
      pixbuf = gdk_pixbuf_new_from_xpm_data((const gchar **) tools_xpm);
      name = "image_tools";
      break;

    default:
      pixbuf = NULL;
      name = NULL;
    }

/* add the entry */
  if (name && pixbuf)
    g_hash_table_insert(sysenv.image_table, name, pixbuf);
  }
}

