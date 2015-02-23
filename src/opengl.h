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

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

/**************/
/* structures */
/**************/

struct canvas_pak
{
gint active;
gint resize;
gpointer model;
gint x,y;
gint width,height;
gint size;
/* openGL transformation data */
GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];
};

struct pipe_pak
{
gdouble v1[3];
gdouble v2[3];
gdouble radius;
gdouble colour[4];
};

/**************/
/* prototypes */
/**************/

struct canvas_pak *gl_new_canvas(gint width, gint height);

void draw_objs(struct canvas_pak *, struct model_pak *);
void gl_draw(struct canvas_pak *, struct model_pak *);
gint gl_init_visual(void);
void gl_init_projection(struct canvas_pak *, struct model_pak *);
void gl_select_box(GtkWidget *);
void gl_project(gdouble *, gint, gint, struct canvas_pak *);
void gl_unproject(gint *, gdouble *, struct canvas_pak *);
gint gl_vertex_visible(gdouble *);
void gl_vertex_window(gint, gint, struct canvas_pak *);

gpointer gl_seek_core(gint, gint, struct model_pak *);
gpointer gl_seek_bond(gint, gint, struct model_pak *);

void set_colour(GdkColor *, gint );
gchar *get_mode_label(struct model_pak *);

void gl_free_points(struct point_pak *);
void gl_init_sphere(struct point_pak *, struct model_pak *);
void gl_draw_sphere(struct point_pak *, gdouble *, gdouble);
void gl_init_circle(struct point_pak *, gint, struct model_pak *);
void gl_draw_circle(struct point_pak *, gdouble *, gdouble);
void gl_draw_ring(struct point_pak *, gdouble *, gdouble, gdouble);
void gl_draw_cylinder(gdouble *, gdouble *, gdouble, guint);
void gl_draw_box(gint, gint, gint, gint, struct canvas_pak *);

void draw_cone(gdouble *, gdouble *, gdouble, gint);
void draw_vector(gdouble *, gdouble *, gdouble);
void draw_arc(gdouble *, gdouble *, gdouble *);

void gl_print_window(gchar *, gint, gint, struct canvas_pak *);
void gl_print_world(gchar *, gdouble, gdouble, gdouble);

void gl_font_free(void);

void stereo_init_window(struct canvas_pak *);
gint stereo_expose_event(GtkWidget *, GdkEventExpose *);
void stereo_close_window(void);
void stereo_open_window(void);
void stereo_draw(void);

void graph_draw(struct canvas_pak *, struct model_pak *);


