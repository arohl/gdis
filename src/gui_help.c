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
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "gui_shorts.h"
#include "dialog.h"
#include "scan.h"
#include "parse.h"
#include "interface.h"
#include "gui_image.h"

extern struct sysenv_pak sysenv;

#define MANUAL_FILE "gdis.manual"

GtkTextTag *help_tag_bold, *help_tag_italic, *help_tag_invisible, *help_tag_fixed;

/**********************/
/* add a single topic */
/**********************/
void help_topic_new(gchar *topic, gchar *text)
{
if (topic)
  {
  if (text)
    g_hash_table_insert(sysenv.manual, g_strdup(topic), g_strdup(text));
  else
    g_hash_table_insert(sysenv.manual, g_strdup(topic), "<empty>\n");
  }
else
  g_hash_table_insert(sysenv.manual, "error", "manual file was probably not found\n");
}

/**************************************************/
/* load the help file and process the manual text */
/**************************************************/
void help_init(void)
{
gint n;
gchar *filename, *line, *topic;
gpointer scan;
GString *buffer;

/* hash table with topic label as the key */
sysenv.manual = g_hash_table_new(g_str_hash, g_str_equal);

/* open the manual file for scanning */
filename = g_build_filename(sysenv.gdis_path, MANUAL_FILE, NULL);
printf("scanning: %s\n", filename);
scan = scan_new(filename);
g_free(filename);
if (!scan)
  {
  help_topic_new(NULL, NULL);
  return;
  }

/* process the manual file */
topic=NULL;
buffer = g_string_new(NULL);
line = scan_get_line(scan);
while (!scan_complete(scan))
  {
  if (g_ascii_strncasecmp(line, "%topic ", 7) == 0)
    {
/* add any old data */
    if (topic && buffer->len)
      {
      help_topic_new(topic, buffer->str);
      buffer = g_string_set_size(buffer, 0);
      g_free(topic);
      }
/* remove trailing <cr> */
    n = strlen(&line[7]);
    topic = g_strndup(&line[7], n-1);
    }
  else
    {
/* append everything that's not a topic to the current buffer */
    if (topic)
      g_string_append(buffer, line);
    }

  line = scan_get_line(scan);
  }

/* done - add the last topic found */
if (topic)
  help_topic_new(topic, buffer->str);

g_string_free(buffer, TRUE);
g_free(topic);

scan_free(scan);
}

/**********************/
/* add a single entry */
/**********************/
void gui_help_populate(gpointer key, gpointer value, gpointer data)
{
gchar *topic = key;
GtkTreeIter iter;
GtkListStore *list = data;

gtk_list_store_append(list, &iter);
gtk_list_store_set(list, &iter, 0, topic, -1);
}

/****************************************/
/* fill out the available manual topics */
/****************************************/
void gui_help_refresh(GtkListStore *list)
{
g_assert(sysenv.manual != NULL);

gtk_list_store_clear(list);
g_hash_table_foreach(sysenv.manual, gui_help_populate, list);
}

/*********************************************/
/* process and display appropriate help text */
/*********************************************/
void gui_help_topic_show(GtkTextBuffer *buffer, gchar *source)
{
gint i, j, n, stop;
#ifdef UNUSED_BUT_SET
gint start;
#endif
gchar *text;
GdkPixbuf *pixbuf;
GtkTextIter a, b;

/* initialize help tags */
if (!help_tag_bold)
  help_tag_bold = gtk_text_buffer_create_tag(buffer,
                                            "bold",
                                            "weight", PANGO_WEIGHT_BOLD,
                                            NULL);

if (!help_tag_italic)
  help_tag_italic = gtk_text_buffer_create_tag(buffer,
                                              "italic",
                                              "style", PANGO_STYLE_ITALIC,
                                              NULL);

/* invisible tag not implemented as of GTK 2.2 */
/*
if (!help_tag_invisible)
  help_tag_invisible = gtk_text_buffer_create_tag(buffer,
                                                 "invisible",
                                                 "invisible",
                                                  TRUE,
                                                  NULL);
*/

/* FIXME - hack due to above
 * ^^^^^ Unfortunately this hack will not work on a non white background
 * _BUG_ --- OVHPA 
if (!help_tag_invisible)
  {
  GdkColor white = {0, 65535, 65535, 65535};
  help_tag_invisible = gtk_text_buffer_create_tag(buffer,
                                                 "invisible",
                                                 "foreground-gdk",
                                                  &white,
                                                  NULL);
  }
*/
if (!help_tag_fixed)
  {
  help_tag_fixed = gtk_text_buffer_create_tag(buffer,
                                              "fixed",
                                              "font", "Fixed 10",
                                              NULL);
  }

/*
tag_bold = gtk_text_buffer_create_tag(buffer, "bold", "style", PANGO_WEIGHT_OBLIQUE, NULL);
*/

gtk_text_buffer_set_text(buffer, source, -1);

n = strlen(source);
#ifdef UNUSED_BUT_SET
start = 0;
#endif
i = 0;
stop = -1;
while (i<n)
  {
/* search for special tags */
  if (source[i] == '\\')
    {
/* process the tag */
    switch (source[i+1])
      {
/* header */
      case 'h':
/* make the tag invisible */
        gtk_text_buffer_get_iter_at_offset(buffer, &a, i);
        gtk_text_buffer_get_iter_at_offset(buffer, &b, i+2);
/*        gtk_text_buffer_apply_tag(buffer, help_tag_invisible, &a, &b);
 * removed due to above trouble with non white BG _BUG_ -- OVHPA */
	gtk_text_buffer_delete(buffer, &a, &b);
	gtk_text_buffer_insert(buffer,&b,"  ",2);
	gtk_text_buffer_get_iter_at_offset(buffer, &b, i+3);

        a = b;
/* apply the bold tag to the rest of the line */
        stop = -1;
        for (j=i+3 ; j<n ; j++)
          {
          if (source[j] == '\n')
            {
            stop = j;
            break;
            }
          }
        if (stop < 0)
          stop = n;
        gtk_text_buffer_get_iter_at_offset(buffer, &b, stop);
        gtk_text_buffer_apply_tag(buffer, help_tag_bold, &a, &b);
        i = stop;
        break;

/* italics */
      case 'i':
/* make the tag invisible */
        gtk_text_buffer_get_iter_at_offset(buffer, &a, i);
        gtk_text_buffer_get_iter_at_offset(buffer, &b, i+2);
/*        gtk_text_buffer_apply_tag(buffer, help_tag_invisible, &a, &b);
 * removed due to above trouble with non white BG _BUG_ -- OVHPA */
	gtk_text_buffer_delete(buffer, &a, &b);
	gtk_text_buffer_insert(buffer,&b,"  ",2);
	gtk_text_buffer_get_iter_at_offset(buffer, &b, i+3);

        a = b;
/* apply the bold tag to the rest of the line */
        stop = -1;
        for (j=i+3 ; j<n ; j++)
          {
          if (source[j] == '\n')
            {
            stop = j;
            break;
            }
          }
        if (stop < 0)
          stop = n;
        gtk_text_buffer_get_iter_at_offset(buffer, &b, stop);
        gtk_text_buffer_apply_tag(buffer, help_tag_italic, &a, &b);
        i = stop;
        break;

/* fixed width font */
      case 'f':
/* make the tag invisible */
        gtk_text_buffer_get_iter_at_offset(buffer, &a, i);
        gtk_text_buffer_get_iter_at_offset(buffer, &b, i+2);
/*        gtk_text_buffer_apply_tag(buffer, help_tag_invisible, &a, &b);
 * removed due to above trouble with non white BG _BUG_ -- OVHPA */
	gtk_text_buffer_delete(buffer, &a, &b);
	gtk_text_buffer_insert(buffer,&b,"  ",2);
	gtk_text_buffer_get_iter_at_offset(buffer, &b, i+3);

        a = b;
/* apply the fixed tag to the rest of the line */
        stop = -1;
        for (j=i+3 ; j<n ; j++)
          {
          if (source[j] == '\n')
            {
            stop = j;
            break;
            }
          }
        if (stop < 0)
          stop = n;
        gtk_text_buffer_get_iter_at_offset(buffer, &b, stop);
        gtk_text_buffer_apply_tag(buffer, help_tag_fixed, &a, &b);
        i = stop;
        break;

/* picture */
      case 'p':
/* get the picture label code */
        stop = -1; 
        for (j=i+3 ; j<n ; j++)
          {
          if (source[j] == ' ' || source[j] == '\n')
            {
            stop = j;
            break;
            }
          }
        if (stop < 0)
          stop = n;

/* hide the tag */
        gtk_text_buffer_get_iter_at_offset(buffer, &a, i);
        gtk_text_buffer_get_iter_at_offset(buffer, &b, stop);
/*        gtk_text_buffer_apply_tag(buffer, help_tag_invisible, &a, &b);
 * removed due to above trouble with non white BG _BUG_ -- OVHPA */
	gtk_text_buffer_delete(buffer, &a, &b);
	gtk_text_buffer_insert(buffer,&b,"                                             ",stop-i);
	gtk_text_buffer_get_iter_at_offset(buffer, &b, i);

	a = b;
/* CURRENT - delete one char to compensate? */
gtk_text_buffer_get_iter_at_offset(buffer, &b, i+1);
gtk_text_buffer_delete(buffer, &a, &b);

/* insert the pixbuf - counts as one character in the text buffer */
        text = g_strndup(&source[i+3], stop-i-3);
        pixbuf = image_table_lookup(text);
        gtk_text_buffer_insert_pixbuf(buffer, &a, pixbuf);

        i = stop;
        break;
      }
    }
  i++;
  }
}

/***************************/
/* topic selection handler */
/***************************/
void gui_help_topic_selected(GtkTreeSelection *selection, gpointer data)
{
gchar *topic, *text;
GtkTreeIter iter;
GtkTreeModel *treemodel;
GtkTextBuffer *buffer;

/* record selection as the active entry */
if (gtk_tree_selection_get_selected(selection, &treemodel, &iter))
  {
  gtk_tree_model_get(treemodel, &iter, 0, &topic, -1);
  text = g_hash_table_lookup(sysenv.manual, topic);
  buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(data));
  gui_help_topic_show(buffer, text);
  }
}

/*******************************/
/* topic list sorting primitve */
/*******************************/
gint help_sort_topics(GtkTreeModel *treemodel, GtkTreeIter *a, GtkTreeIter *b, gpointer data)
{
gchar *str1, *str2;

gtk_tree_model_get(treemodel, a, 0, &str1, -1);
gtk_tree_model_get(treemodel, b, 0, &str2, -1);

return(strcasecmp(str1, str2));
}

/****************************/
/* display a manual browser */
/****************************/
void gui_help_dialog(void)
{
gpointer dialog;
GtkCellRenderer *r;
GtkTreeViewColumn *c;
GtkListStore *list;
GtkTreeSelection *select;
GtkWidget *window, *swin, *hbox, *tree, *view;

/* create dialog */
dialog = dialog_request(MANUAL, "Manual", NULL, NULL, NULL);
if (!dialog)
  return;
window = dialog_window(dialog);
gtk_window_set_default_size(GTK_WINDOW(window), 800, 500);

/* split pane display */
hbox = gtk_hbox_new(FALSE, 0);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), hbox, TRUE, TRUE, 1);

/* left pane - topics browser */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_NEVER, GTK_POLICY_NEVER);
gtk_box_pack_start(GTK_BOX(hbox), swin, FALSE, TRUE, 0);

/* list */
list = gtk_list_store_new(1, G_TYPE_STRING);
tree = gtk_tree_view_new_with_model(GTK_TREE_MODEL(list));
gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(swin), tree);
r = gtk_cell_renderer_text_new();
c = gtk_tree_view_column_new_with_attributes(" ", r, "text", 0, NULL);
gtk_tree_view_append_column(GTK_TREE_VIEW(tree), c);
gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(tree), FALSE);

/* auto sort the topic list */
gtk_tree_sortable_set_default_sort_func(GTK_TREE_SORTABLE(GTK_TREE_MODEL(list)),
                                        help_sort_topics, NULL, NULL); 
gtk_tree_sortable_set_sort_column_id(GTK_TREE_SORTABLE(GTK_TREE_MODEL(list)),
                                     GTK_TREE_SORTABLE_DEFAULT_SORT_COLUMN_ID,
                                     GTK_SORT_ASCENDING);

/* set the width of the topics pane */
gtk_widget_set_size_request(swin, 15*sysenv.gtk_fontsize, -1);

/* right pane - text */
swin = gtk_scrolled_window_new(NULL, NULL);
gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(swin),
                               GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
gtk_box_pack_start(GTK_BOX(hbox), swin, TRUE, TRUE, 0);
view = gtk_text_view_new();
gtk_text_view_set_editable(GTK_TEXT_VIEW(view), FALSE);
gtk_container_add(GTK_CONTAINER(swin), view);

/* configure the text viewing area */
gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(view), GTK_WRAP_WORD);
/* NB: GTK_JUSTIFY_FILL was not supported at the time this was written */
/*
gtk_text_view_set_justification(GTK_TEXT_VIEW(view), GTK_JUSTIFY_FILL);
*/
gtk_text_view_set_left_margin(GTK_TEXT_VIEW(view), PANEL_SPACING);
gtk_text_view_set_right_margin(GTK_TEXT_VIEW(view), PANEL_SPACING);

/* NB: these are associated with the text buffer and must be reallocated each */
/* time a new text buffer is created */
help_tag_bold=NULL;
help_tag_italic=NULL;
help_tag_invisible=NULL;
help_tag_fixed=NULL;

/* topic selection handler */
select = gtk_tree_view_get_selection(GTK_TREE_VIEW(tree));
gtk_tree_selection_set_mode(select, GTK_SELECTION_SINGLE);
g_signal_connect(G_OBJECT(select), "changed",
                 G_CALLBACK(gui_help_topic_selected),
                 view);

/* terminating button */
gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);

/* populate the manual page */
gui_help_refresh(list);

/* expose the dialog */
gtk_widget_show_all(window);
}

/***************/
/* help dialog */
/***************/
/*
      g_string_assign(title,"OpenGL");
      g_string_sprintf(buff,"GL_RENDERER    %s\n", glGetString(GL_RENDERER));
      g_string_sprintfa(buff,"GL_VERSION       %s\n", glGetString(GL_VERSION));
      g_string_sprintfa(buff,"GL_VENDOR        %s\n", glGetString(GL_VENDOR));
*/
