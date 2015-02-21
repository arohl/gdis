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

gpointer mesch_mat_new(gint, gint);
gpointer mesch_vec_new(gint);
void mesch_m_free(gpointer);
void mesch_m_zero(gpointer);
void mesch_v_zero(gpointer);

void mesch_me_set(gpointer, gint, gint, gdouble);
void mesch_me_add(gpointer, gint, gint, gdouble);
void mesch_me_mul(gpointer, gint, gint, gdouble);
void mesch_ve_set(gpointer, gint, gdouble);

gdouble mesch_me_get(gpointer, gint, gint);
gdouble mesch_ve_get(gpointer, gint);

gint mesch_rows_get(gpointer);
gint mesch_cols_get(gpointer);
gint mesch_dim_get(gpointer);

void mesch_sev_compute(gpointer, gpointer, gpointer);

