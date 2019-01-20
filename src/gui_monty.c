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

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "graph.h"
#include "model.h"
#include "parse.h"
#include "scan.h"
#include "task.h"
#include "matrix.h"
#include "surface.h"
#include "spatial.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"
#include "numeric.h"

#define DEBUG_MONTY_GUI 0 

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

void init_monty_settings(struct model_pak * model)
{
  /* for now we set up manually for the mac
  * eventually this should be set by g_find_program_in_path
  */
  g_free(sysenv.monty_path);
  sysenv.monty_path = g_strdup("/Users/menno/src/Monty2/Monty2_GNU_MAC");
  g_free(model->monty.input_cgf);
  model->monty.input_cgf = g_strdup(model->filename);
/* init has moved to model.c */
  return; 
}

/* Unfortunately GTK+ 2.2 and older don't have GtkFileChooserDialog
   so for them a file dialog without GtkFileChooserDialog is given 
   the disadvantage is that this has no easy file filter */
#if GTK_MAJOR_VERSION >= 2 && GTK_MINOR_VERSION >= 4

void select_cgf_file_dialog(GtkButton *button, gpointer data)
{
  GtkWidget *file_chooser;
  GtkFileFilter *cgf_filter;
  
  struct model_pak *model = sysenv.active_model;
 
  file_chooser = gtk_file_chooser_dialog_new("Select a CGF File",
                                             NULL, /* parent window, seems to work without valid pointer */
                                             GTK_FILE_CHOOSER_ACTION_OPEN,
                                             GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                             GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                             NULL);
  
  cgf_filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(cgf_filter, "*.cgf");
  gtk_file_filter_add_pattern(cgf_filter, "*.CGF");
  
  g_object_set(file_chooser, 
               "filter", cgf_filter,
               NULL);
  
  
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    
    g_free(model->monty.input_cgf);
    model->monty.input_cgf = g_strdup(filename);    
    gui_relation_update_widget(&model->monty.input_cgf);
    g_free (filename);
  }
  
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
  
}

void select_surface_file_dialog(GtkButton *button, gpointer data)
{
  GtkWidget *file_chooser;
  GtkFileFilter *surface_filter;
  struct model_pak *model = sysenv.active_model;
  /* todo: file filter */
  file_chooser = gtk_file_chooser_dialog_new("Select a CGF File",
                                             NULL, /* parent window, seems to work without valid pointer */
                                             GTK_FILE_CHOOSER_ACTION_OPEN,
                                             GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                             GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
                                             NULL);
  
  surface_filter = gtk_file_filter_new();
  gtk_file_filter_add_pattern(surface_filter, "*.monty2");  
  
  if (gtk_dialog_run (GTK_DIALOG (file_chooser)) == GTK_RESPONSE_ACCEPT)
  {
    char *filename;
    
    filename = gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));
    
    g_free(model->monty.input_surface);
    model->monty.input_surface = g_strdup(filename);
    gui_relation_update_widget(&model->monty.input_surface);
    
    g_free (filename);
  }
  
  gtk_widget_destroy (GTK_WIDGET(file_chooser));
}

/* gtk+-2 version < 2.4 */
#else 

GtkFileSelection *file_browser;

void select_input_cgf_filename(GtkButton *button, gpointer data)
{
  gchar *filename;
  
  /* get the filename from the file_browser object */
  g_object_get(file_browser, "filename", &filename, NULL);
  
  struct model_pak *model = sysenv.active_model;
  
  /* free and then set the filename */
  g_free(model->monty.input_cgf);
  model->monty.input_cgf = g_strdup(filename);
  
  /* update the widget displaying the input cgf filename */
  gui_relation_update_widget(&model->monty.input_cgf);
  
  g_free(filename);
}

void select_cgf_file_dialog(GtkButton *button, gpointer data)
{
  /* create a new file browser  */
  file_browser = g_object_new(GTK_TYPE_FILE_SELECTION, NULL);
  
  /* set the title of the browser what about a file filter for CGF's?  */
  gtk_window_set_title(GTK_WINDOW(file_browser), "Select a CGF File");
  
  /* connect signals  */
  g_signal_connect(file_browser->ok_button, "clicked", G_CALLBACK(select_input_cgf_filename), NULL);
  g_signal_connect_swapped(file_browser->ok_button, "clicked", 
                           G_CALLBACK(gtk_widget_destroy), file_browser);
  g_signal_connect_swapped(file_browser->cancel_button, "clicked", 
                           G_CALLBACK(gtk_widget_destroy), file_browser);
  
  /* show window  */
  gtk_widget_show(GTK_WIDGET(file_browser));
}

void select_input_surface_filename(GtkButton *button, gpointer data)
{
  gchar *filename;
  
  /* get the filename from the file_browser */
  g_object_get(file_browser, "filename", &filename, NULL);
  
  struct model_pak *model = sysenv.active_model;
  
  /* free and set the filename selected */
  g_free(model->monty.input_surface);
  model->monty.input_surface = g_strdup(filename);
  
  /* update the widget displaying the input cgf filename */
  gui_relation_update_widget(&model->monty.input_surface);
  
  g_free(filename);
}

void select_surface_file_dialog(GtkButton *button, gpointer data)
{
  /* create a new file browser  */
  file_browser = g_object_new(GTK_TYPE_FILE_SELECTION, NULL);
  
  /* set the title of the browser (how about file filter?  */
  gtk_window_set_title(GTK_WINDOW(file_browser), "Select a Surface File");
  
  /* connect signals  */
  g_signal_connect(file_browser->ok_button, "clicked", G_CALLBACK(select_input_surface_filename), NULL);
  g_signal_connect_swapped(file_browser->ok_button, "clicked", 
                           G_CALLBACK(gtk_widget_destroy), file_browser);
  g_signal_connect_swapped(file_browser->cancel_button, "clicked", 
                           G_CALLBACK(gtk_widget_destroy), file_browser);
  
  /* show window  */
  gtk_widget_show(GTK_WIDGET(file_browser));
}
#endif

/* a callback function to set the energy unit to kcal/mol */
void energy_unit_kcal(GtkButton *button, gpointer data)
{
  struct model_pak *model;
  model = sysenv.active_model;
  
  g_free(model->monty.energy_unit);
  model->monty.energy_unit = g_strdup("kcal/mol");    
}

/* a callback function to set the energy unit to kcal/mol */
void energy_unit_kJ(GtkButton *button, gpointer data)
{
  struct model_pak *model;
  model = sysenv.active_model;
  
  g_free(model->monty.energy_unit);
  model->monty.energy_unit = g_strdup("kJ/mol");  
}

void run_crystal_graph(GtkButton *button, gpointer data)
{
  struct model_pak *model = sysenv.active_model;
  calculate_crystal_graph(model);
}

void monty_crystal_graph_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox, *wdgt;
  gint size = 15;
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
  
  wdgt = gui_direct_spin("Cutoff in a-direction", &data->monty.image_x, 1.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Cutoff in b-direction", &data->monty.image_y, 1.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Cutoff in c-direction", &data->monty.image_z, 1.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = GTK_WIDGET(g_object_new(GTK_TYPE_BUTTON,
                                 "label", "Calculate crystal graph", NULL));
  g_signal_connect(GTK_BUTTON(wdgt), 
                   "clicked", G_CALLBACK(run_crystal_graph), NULL);
  gtk_box_pack_start(GTK_BOX(vbox), wdgt, FALSE, FALSE, 0);
}

/*******************************************/
/* following five functions create the     */
/* five pages in the notebook of the monty */
/* settings dialog                        */
/*******************************************/

void monty_input_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox, *table, *frame_hkl, *frame_dirs, *frame_deltamu, *hbox, *h_sep, *wdgt;

  gint size = 15;
  
  /* the top-level vbox, in which all other widgets are packed */
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
  
  
  /* a table with four widgets, two gui_text_entries showing
    the input cgf filename and the input surface filename
    next to them there are two file dialog buttons for selecting the 
    files 
   */
  
  table = g_object_new(GTK_TYPE_TABLE,
                       "n-rows", 2,
                       "n-columns", 2,
                       "homogeneous", FALSE,
                       NULL);
  
  
  gtk_box_pack_start(GTK_BOX(vbox), table, TRUE, TRUE, 0);
  
  /* as the gui_text_entry and gui_stock_button require a box to pack their 
    widgets, we create a hbox for each of the four widgets, and pack the hboxes
    into the table 
   */
  hbox = gtk_hbox_new(FALSE, 0);
  gui_text_entry("Input CGF Filename", &data->monty.input_cgf, TRUE, FALSE, hbox);/*FIX f556b0*/
  gtk_table_attach_defaults(GTK_TABLE(table), hbox, 0,1,0,1);
  
  hbox = gtk_hbox_new(FALSE, 0);
  gui_text_entry("Input Surface Filename", &data->monty.input_surface, TRUE, FALSE, hbox);/*FIX 0be582*/
  gtk_table_attach_defaults(GTK_TABLE(table), hbox, 0,1,1,2);
  
  hbox = gtk_hbox_new(FALSE, 0);
  gui_stock_button(GTK_STOCK_OPEN, select_cgf_file_dialog, NULL, hbox);/*FIX 84d13d*/
  gtk_table_attach_defaults(GTK_TABLE(table), hbox, 1,2,0,1);
  
  hbox = gtk_hbox_new(FALSE, 0);
  gui_stock_button(GTK_STOCK_OPEN, select_surface_file_dialog, NULL, hbox);/*FIX 304a57*/
  gtk_table_attach_defaults(GTK_TABLE(table), hbox, 1,2,1,2);
  
  h_sep = GTK_WIDGET(g_object_new(GTK_TYPE_HSEPARATOR, NULL));
  gtk_widget_set_size_request(h_sep, -1, size);
  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(h_sep), FALSE, FALSE, 0);
  
  /* up next are three gui_text_windows for entering hkl values
     output directories and driving forces. each gui_text_window
     is packed in a frame showing the gui_text_window's purpose
   */
  
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);
  
  frame_hkl = gtk_frame_new("HKL Values");
  gtk_box_pack_start(GTK_BOX(hbox), frame_hkl, TRUE, TRUE, 0);
  wdgt = gui_text_window(&data->monty.hkls, TRUE);
  gtk_container_add(GTK_CONTAINER(frame_hkl), wdgt);
  
  frame_dirs = gtk_frame_new("Output directories");
  gtk_box_pack_start(GTK_BOX(hbox), frame_dirs, TRUE, TRUE, 0);
  wdgt = gui_text_window(&data->monty.output_dirs, TRUE);
  gtk_container_add(GTK_CONTAINER(frame_dirs), wdgt);
  
  frame_deltamu = gtk_frame_new("Driving forces");
  gtk_box_pack_start(GTK_BOX(hbox), frame_deltamu, TRUE, TRUE, 0);
  wdgt = gui_text_window(&data->monty.supersaturations, TRUE);
  gtk_container_add(GTK_CONTAINER(frame_deltamu), wdgt);
  
  /* next we create a radiobutton group containing two radiobuttons
     to switch between kcal/mol and kJ/mol energies in the crystal graph
   */  
  new_radio_group(0, vbox, TT);
  wdgt = add_radio_button("kcal/mol", (gpointer) energy_unit_kcal, data);
//  if (data->monty.energy_unit == "kcal/mol")
  if(find_in_string("kcal/mol",data->monty.energy_unit) != NULL)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(wdgt), TRUE);
  wdgt = add_radio_button("kJ/mol", (gpointer) energy_unit_kJ, data);
//  if (data->monty.energy_unit == "kJ/mol")
  if(find_in_string("kcal/mol",data->monty.energy_unit) != NULL)
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(wdgt), TRUE);
  
  gui_text_entry("Solvation energy", &data->monty.esolv, TRUE, FALSE, vbox);/*FIX 303d31*/
  
  h_sep = GTK_WIDGET(g_object_new(GTK_TYPE_HSEPARATOR, NULL));
  gtk_widget_set_size_request(h_sep, -1, size);
  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(h_sep), FALSE, FALSE, 0);
  
  hbox = gtk_hbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 0);

  
  gui_direct_check("Run 3D nucleation", &data->monty.run_nucleation, NULL, NULL, hbox);
  gui_direct_spin("Initial nucleus size", &data->monty.nucleus_size, 1.0, DBL_MAX, 1.0, NULL, NULL, hbox);/*FIX 126681*/
  //  gui_direct_check("Run diffusion controlled simulation", &data->monty.run_diffusion, NULL, NULL, vbox);

  
}

void monty_output_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox, *wdgt, *h_sep;
  
  gint size = 15;
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

  /* a text entry for entering the output extension */
  wdgt = gui_text_entry("File extension", &data->monty.output_extension, TRUE, FALSE, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  h_sep = GTK_WIDGET(g_object_new(GTK_TYPE_HSEPARATOR, NULL));
  gtk_widget_set_size_request(h_sep, -1, size);
  gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(h_sep), FALSE, FALSE, 0);
  
  /* four output options to write different types of surface file types */
  
  gui_direct_check("Write Monty surface",  &data->monty.write_surface,  NULL, NULL, vbox);
  gui_direct_check("Write XYZ surface",    &data->monty.write_xyz,      NULL, NULL, vbox);
  gui_direct_check("Write Matlab surface", &data->monty.write_matlab,   NULL, NULL, vbox);
  gui_direct_check("Write MSI surface",    &data->monty.write_msi,      NULL, NULL, vbox);
//  gui_direct_check("Write Gaussian Cube file", &data->monty.write_cube, NULL, NULL, vbox);
  
}

void monty_model_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox, *wdgt;
  gint size = 15;
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
  
  /* a checkbox for selecting the spiral growth mechanism */
  gui_direct_check("Spiral growth", &data->monty.spiral, NULL, NULL, vbox);
  
  /* four spinners selecting the number of x-steps, y-steps, the temperature and the kinetics factor */
  
  wdgt = gui_direct_spin("X-steps", &data->monty.xsteps, -DBL_MAX, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Y-steps", &data->monty.ysteps, -DBL_MAX, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Temperature", &data->monty.temperature, 0.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Kinetics", &data->monty.xsteps, 0.0, 1.0, 0.01, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
}

void monty_monitor_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox;
  
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

  /* three different surface monitoring options */
  gui_direct_check("Create a multi-frame .xyz file", &data->monty.multi_frame_xyz, NULL, NULL, vbox);  
  gui_direct_check("Monitor average height",     &data->monty.monitor_height, NULL, NULL, vbox);
  gui_direct_check("Monitor surface energy",     &data->monty.monitor_energy, NULL, NULL, vbox);
  gui_direct_check("Monitor height correlation", &data->monty.monitor_hhcorr, NULL, NULL, vbox);
  gui_direct_check("Monitor diffusion profile",  &data->monty.monitor_diffusion_profile, NULL, NULL, vbox);

}

void monty_run_box(GtkWidget *box, struct model_pak *data)
{
  GtkWidget *vbox, *wdgt;
  
  /* general width of entry widgets  */
  gint size = 15; 
  
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_box_pack_start(GTK_BOX(box), vbox, TRUE, TRUE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);

  wdgt = gui_text_entry("Random Seed",   &data->monty.random_seed, TRUE, FALSE, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Rows",       &data->monty.rows,      5.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);

  wdgt = gui_direct_spin("Columns",    &data->monty.cols,      5.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Layers",     &data->monty.layers,    5.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Increment",  &data->monty.increment, 2.0, DBL_MAX, 1.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Relaxation", &data->monty.relax,     0.0, DBL_MAX, 1000.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Cycles",     &data->monty.cycles,    0.0, DBL_MAX, 1000.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
  
  wdgt = gui_direct_spin("Moves",      &data->monty.moves,     0.0, DBL_MAX, 1000.0, NULL, NULL, vbox);
  gtk_widget_set_size_request(wdgt, sysenv.gtk_fontsize*size, -1);
}

void write_input_file(gchar *input_file, gpointer data)
{
  FILE *fp;
  struct model_pak *model = data;
  gchar *last_dirsep, *yes, *no;
  gint i;
  
  yes = g_strdup("yes");
  no  = g_strdup("no");

  unlink(input_file );
  
  fp = fopen(input_file, "wt");
  if (!fp)
  {
    gui_text_show(ERROR, "Unable to open Monty input file for writing");
    return;
  }
  
  fprintf(fp, "#Monty input file generated by GDIS %f\n", VERSION);
  fprintf(fp, "&INPUT\n");

  /* remove any whitespace, this also includes trailing \n characters */
  g_strstrip(model->monty.hkls);
  
  /* split the string by endline characters */
  gchar **hkls = g_strsplit(model->monty.hkls, "\n", INT_MAX);  
  
  /* print all hkl values to the file */
  for (i = 0; *(hkls+i); ++i)
  {
    fprintf(fp, " -hkl %s\n", *(hkls+i));
  } 
  
  /* remove any whitespace, this also includes trailing \n characters */
  g_strstrip(model->monty.supersaturations);
  gchar **supersats = g_strsplit(model->monty.supersaturations, "\n", INT_MAX);
  for (i = 0; *(supersats+i); ++i)
  {
    fprintf(fp, " -supersat %s\n", *(supersats +i));
  }
  g_strfreev(supersats);
  
  last_dirsep = g_strrstr(model->monty.input_surface, DIR_SEP);
  
  if (last_dirsep)
    fprintf(fp, " -surface %s\n", (last_dirsep + 1));
  else
  {
    /* TODO: make this work when no complete path is selected  */
    fprintf(fp, "# -surface\n");
  }
  
  last_dirsep = g_strrstr(model->monty.input_cgf, DIR_SEP);
  
  fprintf(fp, " -cgf %s\n", (last_dirsep + 1));
  fprintf(fp, " -directory %s\n", g_strndup(model->monty.input_cgf, (last_dirsep - model->monty.input_cgf)));
  fprintf(fp, " -energy_unit %s\n", model->monty.energy_unit);
  fprintf(fp, " -esolv %s\n", model->monty.esolv);
  
  fprintf(fp, "&END\n\n&OUTPUT\n");
  
  /* remove any whitespace, this also includes trailing \n characters */
  g_strstrip(model->monty.output_dirs);
  gchar **out_dirs = g_strsplit(model->monty.output_dirs, "\n", INT_MAX);
  for (i = 0; *(out_dirs + i); ++i)
  {
    fprintf(fp, " -directory %s %s\n", *(hkls + i), *(out_dirs + i));
  }
  
  g_strfreev(hkls);
  g_strfreev(out_dirs);
  
  fprintf(fp, " -extension %s\n", model->monty.output_extension);
  fprintf(fp, " -surface %s\n", (model->monty.write_surface ? yes : no));
  fprintf(fp, " -xmm %s\n", (model->monty.write_xyz ? yes : no));
  fprintf(fp, " -xyz %s\n", (model->monty.write_matlab ? yes : no));
  fprintf(fp, " -msi %s\n", (model->monty.write_msi ? yes : no));
  
  fprintf(fp, "&END\n\n&MODEL\n");
  
  fprintf(fp, " -spiral %s\n", (model->monty.write_matlab ? yes : no));
  fprintf(fp, " -xsteps %d\n", nearest_int(model->monty.xsteps));
  fprintf(fp, " -ysteps %d\n", nearest_int(model->monty.ysteps));
  fprintf(fp, " -kinetics %2.2f\n", model->monty.kinetics);
  fprintf(fp, " -temperature %10.2f\n", model->monty.temperature);

  fprintf(fp, "&END\n\n&MONITOR\n");
  
  fprintf(fp, " -height %s\n", (model->monty.monitor_height ? yes : no));
  fprintf(fp, " -energy %s\n", (model->monty.monitor_energy ? yes : no));
  fprintf(fp, " -hhcorr %s\n", (model->monty.monitor_hhcorr ? yes : no));
  fprintf(fp, " -gdis %s\n", (model->monty.multi_frame_xyz ? yes : no));

  fprintf(fp, "&END\n\n&RUN\n");
  
  fprintf(fp, " -randomseed %s\n", model->monty.random_seed);
  fprintf(fp, " -rows %d\n", nearest_int(model->monty.rows));
  fprintf(fp, " -cols %d\n", nearest_int(model->monty.cols));
  fprintf(fp, " -layers %d\n", nearest_int(model->monty.layers));
  fprintf(fp, " -increment %d\n", nearest_int(model->monty.increment));
  fprintf(fp, " -relaxation %d\n", nearest_int(model->monty.relax));
  fprintf(fp, " -cycles %d\n", nearest_int(model->monty.cycles));
  fprintf(fp, " -moves %d\n", nearest_int(model->monty.moves));
  
  fprintf(fp, "&END\n");
  g_free(yes);
  g_free(no);
  fclose(fp);
}

gint exec_monty(const gchar *input)
{
  gint status=0;
  gchar *cmd;
  
  /* checks */
  if (!sysenv.monty_path)
    return(-1);
  
  /* delete the old file to be sure output is only data from current run */
#if _WIN32
  /* NB: must enclose paths with single quotes */
  cmd = g_strdup_printf("'%s'", sysenv.monty_path);
#else
  cmd = g_strdup_printf("%s  %s", sysenv.monty_path, input);
#endif
  
#if DEBUG_MONTY_GUI
  printf("executing: [%s]\n",cmd);
#endif
  
  task_sync(cmd);
  
  /* done */
  g_free(cmd);
  return(status);
}

void exec_monty_task(gpointer ptr, gpointer data)
{
    gchar *inpfile;
    struct model_pak *model = ptr;
    struct task_pak *task = data;
    
    /* checks */
    g_assert(model != NULL);
    g_assert(task != NULL);
    
    /* construct fullpath input filename - required for writing */
#if _WIN32
    inpfile = g_build_filename(sysenv.cwd, "input.monty2", NULL);
    /* no status file for the moment, as there are multiple output files for a
      monty simulation, so unable to follow a single one *
    task->status_file = g_build_filename(sysenv.cwd, "gulp.got", NULL);
    */
#else
    inpfile = g_build_filename(sysenv.cwd, "input.monty2", NULL);
    /* no status file for the moment, as there are multiple output files for a
      monty simulation, so unable to follow a single one *
    task->status_file = g_build_filename(sysenv.cwd, model->gulp.out_file, NULL);
     */
#endif
    
#if DEBUG_MONTY_GUI
    printf(" input file: %s\n", inpfile);

#endif
    
    write_input_file(inpfile, model);
    g_free(inpfile);
    
    exec_monty("input.monty2");
}

void gui_monty_task(GtkWidget *w, gpointer data)
{
  if (!sysenv.monty_path)
  {
    gui_text_show(ERROR, "Monty executable was not found.\n");
    return;
  }
  /* put a task that runs monty. for now, no postprocessing */
  task_new("Monty", &exec_monty_task, data, NULL, NULL, data);
}

void monty_execute(GtkWidget *w, gpointer data)
{
  struct model_pak *model = sysenv.active_model;
  if (model) 
  {       
#if DEBUG_MONTY_GUI
    
    /* &INPUT */
    g_message("HKLs %s", model->monty.hkls);
    g_message("Driving forces %s", model->monty.supersaturations);
    g_message("Input surface %s", model->monty.input_surface);
    /* TODO: use DIRSEP */
    g_message("Input directory %s", model->monty.input_dir);
    g_message("Input CGF %s", model->monty.input_cgf);
    g_message("Energy unit %s", model->monty.energy_unit);
    g_message("Esolv %s", model->monty.esolv);
    
    /* &OUTPUT */
    g_message("Output directories %s", model->monty.output_dirs);
    g_message("Output extension %s", model->monty.output_extension);
    g_message("Write surface %s", model->monty.write_surface ? "yes" : "no");
    g_message("Write xyz %s", model->monty.write_xyz ? "yes" : "no");
    g_message("Write matlab %s",model->monty.write_matlab ? "yes" : "no");
    g_message("Write MSI %s", model->monty.write_msi ? "yes" : "no");
    
    /* &MODEL */
    g_message("Spiral growth %s", model->monty.spiral ? "yes" : "no");
    g_message("x steps %2f", model->monty.xsteps);
    g_message("y steps %2f", model->monty.ysteps);
    g_message("kinetics %2f", model->monty.kinetics);
    g_message("temperature %2f", model->monty.temperature);
    
    /* &MONITOR */
    g_message("Monitor height %s", model->monty.monitor_height ? "yes" : "no");
    g_message("Monitor energy %s", model->monty.monitor_energy ? "yes" : "no");
    g_message("Monitor hhcorr %s", model->monty.monitor_hhcorr ? "yes" : "no");  
    
    /* &RUN */
    g_message("random seed %s", model->monty.random_seed);
    g_message("Number rows %10f", model->monty.rows); 
    g_message("Number cols %10f", model->monty.cols);
    g_message("Number layers %10f", model->monty.layers);
    g_message("Increment %10f", model->monty.increment);
    g_message("Relaxations %10f", model->monty.relax);
    g_message("Cycles %10f", model->monty.cycles);
    g_message("Moves %10f", model->monty.moves);
#endif
    gui_monty_task(NULL, model);
  }
  else 
    g_message("no model active");
}

void gui_monty_widget(GtkWidget *box, struct model_pak *data)
{
  GString *line;
  GtkWidget *page;
  GtkWidget *label, *notebook;
  
  /* string manipulation scratchpad */
  line = g_string_new(NULL);
  
  /* create notebook */
  notebook = gtk_notebook_new();
  gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);
  gtk_container_add(GTK_CONTAINER(box), notebook);
  gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), FALSE);  
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Crystal Graph ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_crystal_graph_box(page, data);
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Input ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_input_box(page, data);
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Output ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_output_box(page, data);
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Model ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_model_box(page, data);
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Monitor ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_monitor_box(page, data);
  
  page = gtk_vbox_new(FALSE,0);
  label = gtk_label_new (" Run ");
  gtk_notebook_append_page(GTK_NOTEBOOK(notebook), page, label);
  monty_run_box(page, data);
  
  /* done */
  gtk_widget_show_all(box);
  
  g_string_free(line, TRUE);
  gui_model_select(data);
  
}

void monty_dialog(void)
{
  gpointer dialog;
  GtkWidget *window, *frame, *vbox;
  struct model_pak *model;
  
  model = sysenv.active_model;
  if (!model)
    return;
  if (!model->monty.input_cgf)
    init_monty_settings(model);
  
  /* request a MONTY dialog */
  dialog = dialog_request(MONTY, "MONTY configuration", NULL, NULL, model);
  if (!dialog)
    return;
  window = dialog_window(dialog);
  
  frame = gtk_frame_new(NULL);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, FALSE, FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
  
  vbox = gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
  gui_monty_widget(vbox, model);
  
  /* terminating button */
  gui_stock_button(GTK_STOCK_EXECUTE, monty_execute, model,
                     GTK_DIALOG(window)->action_area);
  
  gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                     GTK_DIALOG(window)->action_area);
  
  gtk_widget_show_all(window);
}
