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

/*************/
/* constants */
/*************/

enum {HF_CURVATURE, HF_LAST};
enum {SPATIAL_LINE, SPATIAL_SURFACE, SPATIAL_SOLID};	
enum {SPATIAL_GENERIC, SPATIAL_VECTOR, SPATIAL_MORPHOLOGY};	

/**************/
/* structures */
/**************/

struct spatial_pak
{
gint type;             /* SPATIAL_GENERIC, or special type */
gint size;             /* vertices per spatial method primitive */
gint method;           /* openGL drawing */
gint material;         /* 1 = wire, 2 = surface (wire or solid), 3 = solid */
gint periodic;         /* create images */
gint show_label;
gchar *label;
gdouble x[3];          /* label position */
gdouble c[3];          /* label colour */
gpointer data;         /* general purpose associated data */
GSList *list;          /* vertex list */
};

/**************/
/* prototypes */
/**************/

void spatial_destroy(gpointer, struct model_pak *);
void spatial_destroy_by_type(gint, struct model_pak *);
void spatial_destroy_by_label(const gchar *, struct model_pak *);
void spatial_destroy_all(struct model_pak *);

void delete_vector_at(gdouble *, struct model_pak *);

gpointer spatial_new(const gchar *, gint , gint, gint, struct model_pak *);
void spatial_vertex_add(gdouble *, gdouble *, gpointer);
void spatial_vnorm_add(gdouble *, gdouble *, gdouble *, gpointer);
void spatial_point_add(struct core_pak *, struct model_pak *);
gpointer spatial_build_facet(gdouble *, GSList *, struct model_pak *);

void spatial_vdata_add(gdouble *, gdouble *, gdouble *, gpointer, gpointer);

void spatial_data_set(gpointer, gpointer);
gpointer spatial_data_get(gpointer);

void compute_plane(gdouble *, GList *); 
void compute_normal(gdouble *, gdouble *, gdouble *, gdouble *); 

