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

/* text manipulation prototypes */
gint copy_items(gchar *, gchar *, gint, gint);
void capitals(gchar *, gchar *);
gchar **get_tokens(gchar *, gint);
gchar *get_token_pos(gchar *, gint);
gchar **tokenize(const gchar *, gint *);
gchar *find_char(const gchar *, gint, gint);
gint char_count(const gchar *, gchar);
gdouble str_to_float(const gchar *);
gint str_is_float(const gchar *);
gint *get_keyword(gchar *, gint);
gint num_keys(void);
gint char_count(const gchar *, gchar);
GSList *get_keywords(gchar *);
GSList *get_keywords_anywhere(gchar *);
void strip_extra(gchar *);
gchar *parse_extension_set(const gchar *, const gchar *);
void parse_space_replace(gchar *, gchar);
gchar *parse_strip_extension(const gchar *);
gchar *parse_strip_newline(const gchar *);
gchar *parse_strip(const gchar *);

gchar *parse_getline_hidden(void);

