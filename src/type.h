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

gpointer type_new(void);

void type_print(gpointer);
void type_free(gpointer);
void type_ff_set(gint, const gchar *, gpointer);
void type_charge_set(gint, gdouble, gpointer);
void type_rule_add(gint, gint, const gchar *, gpointer);
void type_dreiding_gasteiger(struct model_pak *, gint);

gint type_check_list(GSList *);
gint type_apply(gpointer, GSList *);

