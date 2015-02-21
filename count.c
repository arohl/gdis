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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gdis.h"
#include "count.h"
#include "numeric.h"

/**************************/
/* construction primitive */
/**************************/
#define DEBUG_COUNT_NEW 0
gpointer count_new(gdouble start, gdouble stop, gdouble step)
{
struct count_pak *count;

g_assert(stop > step);
g_assert(stop-step >= step);
g_assert(step >= 0.0);

count = g_malloc(sizeof(struct count_pak));

/* copy extents */
count->start = start;
count->stop = stop;
count->step = step;

/* compute and allocate bin array */
count->size = nearest_int (0.5 + ((stop-start) / step) );
count->bins = g_malloc0(count->size * sizeof(gint));

#if DEBUG_COUNT_NEW
printf("new count from: %f - %f, step: %f, size: %d\n", start, stop, step, count->size); 
#endif

return(count);
}

/*************************/
/* destruction primitive */
/*************************/
void count_free(gpointer data)
{
struct count_pak *count = data;

g_assert(count != NULL);

g_free(count->bins);
g_free(count);
}

/***************************/
/* debug/contents function */
/***************************/
void count_stats(gpointer data)
{
gint i, value, sum, nonzero;
struct count_pak *count = data;

g_assert(count != NULL);

sum = nonzero = 0;
for (i=count->size ; i-- ; )
  {
  value = *(count->bins+i);
  sum += value;
  }

printf("Count %p: size=%d : sum=%d\n", count, count->size, sum);
}

/*************************/
/* data access primitive */
/*************************/
gint *count_bins(gpointer data)
{
struct count_pak *count = data;

return(count->bins);
}

gint count_size(gpointer data)
{
struct count_pak *count = data;

return(count->size);
}

gdouble count_stop(gpointer data)
{
struct count_pak *count = data;

return(count->stop);
}

/**************************/
/* attempt to bin a value */
/**************************/
/* NB: returns -1, 0, 1 if value is less, within, or outside the count range */
gint count_insert(gdouble value, gpointer data)
{
gint n;
gdouble offset;
struct count_pak *count = data;

g_assert(count != NULL);

if (value < count->start)
  return(-1);
if (value > count->stop)
  return(1);

offset = value - count->start;
offset /= count->step;

n = nearest_int(offset);

g_assert(n >= 0);

if (n >= count->size)
  {
printf("BAD VALUE = %f : size = %d\n", value, n);
  }

g_assert(n < count->size);

*(count->bins+n) += 1; 

return(0);
}

/*********************************************************/
/* assuming congruent, merge the bin data for two counts */
/*********************************************************/
void count_add(gpointer a, gpointer b)
{
gint i;
struct count_pak *c1=a, *c2=b;

g_assert(c1 != NULL);
g_assert(c2 != NULL);

/* TODO - can get more fancy with the merge st mismatch range / size etc are handled */
g_assert(c1->size == c2->size);

for (i=c1->size ; i-- ; )
  *(c1->bins+i) += *(c2->bins+i);
}

