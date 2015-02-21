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

struct model_pak *model_new(void);
void model_init(struct model_pak *);
gint model_prep(struct model_pak *);
void model_free(struct model_pak *);
void model_delete(struct model_pak *);
gpointer model_dup(struct model_pak *);

void model_content_refresh(struct model_pak *);

void property_free(gpointer);
void property_add_ranked(guint, 
                         const gchar *,
                         const gchar *,
                         struct model_pak *);
guint property_rank(gpointer);
gchar *property_label(gpointer);
gchar *property_value(gpointer);
gchar *property_lookup(gchar *, struct model_pak *);

void gulp_init(gpointer);
void gamess_init(gpointer);

