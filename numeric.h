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

/* prototypes */

void init_math(void);

gdouble tbl_sin(gdouble);
gdouble tbl_cos(gdouble);

void ca_deg(gdouble *);
void ca_rad(gdouble *);

gdouble angle_x_compute(gdouble, gdouble);

gdouble decimal_round(gdouble, gint);
gint nearest_int(gdouble);

gint gcd(gint, gint);

gdouble erf(gdouble);
gdouble erfc(gdouble);

gulong mytimer(void);

void sort(gint, gdouble *);

gdouble min(gint, gdouble *);
gdouble max(gint, gdouble *);

void spline(double *, double *, int, double, double, double *);
void splint(double *, double *, double *, int, double, double *);
void fft(gdouble *, gint, gint);

