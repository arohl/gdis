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

/********/
/* data */
/********/

struct defect_pak
{
/* debugg option */
gint cluster;
/* Mott-Littleton */
gint cleave;
gint neutral;
gdouble region[2];
/* geometry */
gdouble orient[3];
gdouble burgers[3];
gdouble origin[2];
gdouble center[2];
};

/**************/
/* prototypes */
/**************/

void defect_new(struct defect_pak *, struct model_pak *);

