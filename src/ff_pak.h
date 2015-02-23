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

#define FF_MAX_ATOMS 4
#define FF_MAX_SYMBOL 8
#define FF_MAX_DATA 7

struct forcefield_pak
{
/* FF type */
gint type;
/* atomic symbols */
gchar atom[FF_MAX_ATOMS][FF_MAX_SYMBOL];

/* units for forcefield parameters */
gint bond_units;
gint data_units;

/* FF parameters */
gdouble bond_value;
gdouble data[FF_MAX_DATA];

/* internal parsing data */
gint atoms_expected;
gint atoms_current;
gint data_expected;
gint data_current;
gint bond_expected;

/* indicates which of the data[] values should contain the bond length/angle (if any) */
gint bond_index;
};

/* prototypes */
gint ff_match_code(struct forcefield_pak *, gint *, gint);
gint ff_match_label(struct forcefield_pak *, gchar **, gint);

void ff_swap_atoms(struct forcefield_pak *, gint, gint);

