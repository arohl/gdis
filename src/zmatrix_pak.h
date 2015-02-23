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

struct zmat_pak
{
/* TODO - have 2 lists - the 1st includes the zval / variable format */
/* the second is the actual core_pak list subset - for easy updating */
gint fractional;
gdouble distance_scale;
gdouble angle_scale;
gchar *distance_units;
gchar *angle_units;
GSList *zlines;
GSList *zcores;
GHashTable *vars;
GHashTable *consts;
};

struct zval_pak
{
gint type;
gchar *elem;
gint connect[3];
gint fitting[3];
gpointer name[3];   /* if NULL -> value only ie not a variable */
gdouble value[3];
};

