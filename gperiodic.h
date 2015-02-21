/* 
   Copyright (C) 1999 Kyle R. Burton, All Rights Reserved
   mortis@voicenet.com
   http://www.bgw.org
   http://www.voicenet.com/~mortis

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

/* NEW - simplfied struct (all element data is got from the unified database) */

struct table_entry {
  int   x;                /* x position in the table */
  int   y;                /* y position in the table */
  GtkStyle    *style;
  GtkWidget   *button;
  GtkTooltips *tooltip;
};

struct table_entry table[] = {
{ 1,  1, NULL, NULL, NULL},
{ 18,  1, NULL, NULL, NULL},
{ 1,  2, NULL, NULL, NULL},
{ 2,  2, NULL, NULL, NULL},
{ 13,  2, NULL, NULL, NULL},
{ 14,  2, NULL, NULL, NULL},
{ 15,  2, NULL, NULL, NULL},
{ 16,  2, NULL, NULL, NULL},
{ 17,  2, NULL, NULL, NULL},
{ 18,  2, NULL, NULL, NULL},
{ 1,  3, NULL, NULL, NULL},
{ 2,  3, NULL, NULL, NULL},
{ 13,  3, NULL, NULL, NULL},
{ 14,  3, NULL, NULL, NULL},
{ 15,  3, NULL, NULL, NULL},
{ 16,  3, NULL, NULL, NULL},
{ 17,  3, NULL, NULL, NULL},
{ 18,  3, NULL, NULL, NULL},
{ 1,  4, NULL, NULL, NULL},
{ 2,  4, NULL, NULL, NULL},
{ 3,  4, NULL, NULL, NULL},
{ 4,  4, NULL, NULL, NULL},
{ 5,  4, NULL, NULL, NULL},
{ 6,  4, NULL, NULL, NULL},
{ 7,  4, NULL, NULL, NULL},
{ 8,  4, NULL, NULL, NULL},
{ 9,  4, NULL, NULL, NULL},
{ 10,  4, NULL, NULL, NULL},
{ 11,  4, NULL, NULL, NULL},
{ 12,  4, NULL, NULL, NULL},
{ 13,  4, NULL, NULL, NULL},
{ 14,  4, NULL, NULL, NULL},
{ 15,  4, NULL, NULL, NULL},
{ 16,  4, NULL, NULL, NULL},
{ 17,  4, NULL, NULL, NULL},
{ 18,  4, NULL, NULL, NULL},
{ 1,  5, NULL, NULL, NULL},
{ 2,  5, NULL, NULL, NULL},
{ 3,  5, NULL, NULL, NULL},
{ 4,  5, NULL, NULL, NULL},
{ 5,  5, NULL, NULL, NULL},
{ 6,  5, NULL, NULL, NULL},
{ 7,  5, NULL, NULL, NULL},
{ 8,  5, NULL, NULL, NULL},
{ 9,  5, NULL, NULL, NULL},
{ 10,  5, NULL, NULL, NULL},
{ 11,  5, NULL, NULL, NULL},
{ 12,  5, NULL, NULL, NULL},
{ 13,  5, NULL, NULL, NULL},
{ 14,  5, NULL, NULL, NULL},
{ 15,  5, NULL, NULL, NULL},
{ 16,  5, NULL, NULL, NULL},
{ 17,  5, NULL, NULL, NULL},
{ 18,  5, NULL, NULL, NULL},
{ 1,  6, NULL, NULL, NULL},
{ 2,  6, NULL, NULL, NULL},
{ 3,  6, NULL, NULL, NULL},
{ 4,  9, NULL, NULL, NULL},
{ 5,  9, NULL, NULL, NULL},
{ 6,  9, NULL, NULL, NULL},
{ 7,  9, NULL, NULL, NULL},
{ 8,  9, NULL, NULL, NULL},
{ 9,  9, NULL, NULL, NULL},
{ 10,  9, NULL, NULL, NULL},
{ 11,  9, NULL, NULL, NULL},
{ 12,  9, NULL, NULL, NULL},
{ 13,  9, NULL, NULL, NULL},
{ 14,  9, NULL, NULL, NULL},
{ 15,  9, NULL, NULL, NULL},
{ 16,  9, NULL, NULL, NULL},
{ 17,  9, NULL, NULL, NULL},
{ 4,  6, NULL, NULL, NULL},
{ 5,  6, NULL, NULL, NULL},
{ 6,  6, NULL, NULL, NULL},
{ 7,  6, NULL, NULL, NULL},
{ 8,  6, NULL, NULL, NULL},
{ 9,  6, NULL, NULL, NULL},
{ 10,  6, NULL, NULL, NULL},
{ 11,  6, NULL, NULL, NULL},
{ 12,  6, NULL, NULL, NULL},
{ 13,  6, NULL, NULL, NULL},
{ 14,  6, NULL, NULL, NULL},
{ 15,  6, NULL, NULL, NULL},
{ 16,  6, NULL, NULL, NULL},
{ 17,  6, NULL, NULL, NULL},
{ 18,  6, NULL, NULL, NULL},
{ 1,  7, NULL, NULL, NULL},
{ 2,  7, NULL, NULL, NULL},
{ 3,  7, NULL, NULL, NULL},
{ 4,  10, NULL, NULL, NULL},
{ 5,  10, NULL, NULL, NULL},
{ 6,  10, NULL, NULL, NULL},
{ 7,  10, NULL, NULL, NULL},
{ 8,  10, NULL, NULL, NULL},
{ 9,  10, NULL, NULL, NULL},
{ 10,  10, NULL, NULL, NULL},
{ 11,  10, NULL, NULL, NULL},
{ 12,  10, NULL, NULL, NULL},
{ 13,  10, NULL, NULL, NULL},
{ 14,  10, NULL, NULL, NULL},
{ 15,  10, NULL, NULL, NULL},
{ 16,  10, NULL, NULL, NULL},
{ 17,  10, NULL, NULL, NULL},
{ 4,  7, NULL, NULL, NULL},
{ 5,  7, NULL, NULL, NULL},
{ 6,  7, NULL, NULL, NULL},
{ 7,  7, NULL, NULL, NULL},
{ 8,  7, NULL, NULL, NULL},
{ 9,  7, NULL, NULL, NULL},
{ 10,  7, NULL, NULL, NULL},
{ 11,  7, NULL, NULL, NULL},
{ 12,  7, NULL, NULL, NULL},
{ 13,  7, NULL, NULL, NULL},
{ 14,  7, NULL, NULL, NULL},
{ 15,  7, NULL, NULL, NULL},
{ 16,  7, NULL, NULL, NULL},
{ 17,  7, NULL, NULL, NULL},
{ 18,  7, NULL, NULL, NULL}
};


void refresh_table(void);

