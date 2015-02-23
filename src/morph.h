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

enum {DHKL, EQUIL_UN, EQUIL_RE, GROWTH_UN, GROWTH_RE, MORPH_BBPA};

/**************/
/* structures */
/**************/

struct vertex_pak
{
/* coords */
gdouble x[3];
gdouble rx[3];
gdouble n[3];
/* adjacencies */
GSList *adj;
};

/**************/
/* prototypes */
/**************/

gint morph_build(struct model_pak *);
gint morph_create(gchar *);
gint facet_visible(struct model_pak *, struct plane_pak *);
gint facet_equiv(struct model_pak *, gint *, gint *);

gint calc_valid_shifts(struct model_pak *, struct plane_pak *);

