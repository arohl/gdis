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

/******************/
/* MDI structures */
/******************/
/* TODO - combine box_pak & cand_pak? */
struct box_pak
{
gint component;
gint x;
gint y;
gint z;
};

struct mdi_pak
{
gint box_dim;
gdouble latt_sep;
gint num_comp;
gint *comp_idx;
gint *comp_req;
gint *comp_done;
gint *array;
};

struct cand_pak
{
gint pos;
gint x;
gint y;
gint z;
};

