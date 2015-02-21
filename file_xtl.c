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
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "interface.h"

#define DEBUG_MORE 0
#define MAX_KEYS 15

/* Added by C. Fisher 2004 */
#define CELL_LENGTH 10.0 /* Box length for non-periodic system */

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/**********************/
/* save in XTL format */
/**********************/
gint write_xtl(gchar *filename, struct model_pak *data)
{
gint i;
gint invert=FALSE;
gint charge;
gint num_tokens;
gchar **buff, *t, *symbol;
gdouble vec[3], tmat[9];
GSList *list;
struct elem_pak elem_data;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
/* g_return_val_if_fail(data->periodic, 1); */

/* open the file */
fp = fopen(filename,"w");
if (!fp)
  return(1);

gui_text_show(STANDARD, "Saving file in XTL format!\n");

/* header */
if( data->title)
  fprintf(fp,"TITLE %s\n",data->title);
if( data->periodic > 0 )
  {
  fprintf(fp,"DIMENSION %d\n",data->periodic);
  fprintf(fp,"CELL\n");
  if (data->periodic == 3)
    {
    fprintf(fp,"%f  %f  %f  %f  %f  %f\n",data->pbc[0],data->pbc[1],data->pbc[2],
          R2D*data->pbc[3], R2D*data->pbc[4], R2D*data->pbc[5]);

    fprintf(fp,"SYMMETRY  NUMBER %d", data->sginfo.spacenum);
    /* Write space group after removing spaces */
    if (data->sginfo.spacename)
      {
      buff = tokenize(data->sginfo.spacename, &num_tokens);
      fprintf(fp,"  LABEL ");
      for( i=0; i < num_tokens; i++)
        {
        t = *(buff+i);
        while( *t)
          {
          fprintf(fp,"%c",toupper(*t));
          t++;
          }
        }
      }

     if( data->sginfo.spacenum > 1 && data->sginfo.spacenum < 16 )
       {
       if( data->pbc[4] != 90.0/R2D )
         fprintf(fp,"  QUALIFIER B_UNIQUE");
       else if( data->pbc[3] != 90.0/R2D )
         fprintf(fp,"  QUALIFIER A_UNIQUE");
       else if( data->pbc[5] != 90.0/R2D )
         fprintf(fp,"  QUALIFIER C_UNIQUE");
       }

    if (data->sginfo.cellchoice)
      {
      if( data->sginfo.spacenum == 146 || data->sginfo.spacenum == 148 ||
          data->sginfo.spacenum == 155 ||
          data->sginfo.spacenum == 160 || data->sginfo.spacenum == 161 ||
          data->sginfo.spacenum == 166 || data->sginfo.spacenum == 167 )
            fprintf(fp,"  QUALIFIER %s", data->sginfo.cellchoice==1?"HEXAGONAL":"RHOMBOHEDRAL");
      else if( data->sginfo.spacenum == 48 || data->sginfo.spacenum == 70 ||
          data->sginfo.spacenum == 85 || data->sginfo.spacenum == 86 ||
          data->sginfo.spacenum == 88 || data->sginfo.spacenum == 125 ||
          data->sginfo.spacenum == 126 || data->sginfo.spacenum == 129 ||
          data->sginfo.spacenum == 130 || data->sginfo.spacenum == 133 ||
          data->sginfo.spacenum == 134 || data->sginfo.spacenum == 137 ||
          data->sginfo.spacenum == 138 || data->sginfo.spacenum == 141 ||
          data->sginfo.spacenum == 142 || data->sginfo.spacenum == 201 ||
          data->sginfo.spacenum == 203 || data->sginfo.spacenum == 222 ||
          data->sginfo.spacenum == 224 || data->sginfo.spacenum == 227 ||
          data->sginfo.spacenum == 228 )
            fprintf(fp,"  QUALIFIER ORIGIN_%d", data->sginfo.cellchoice);
      else if( data->sginfo.spacenum == 50 || data->sginfo.spacenum == 59 ||
                data->sginfo.spacenum == 68 )
           {
           if( data->sginfo.cellchoice % 2 == 1)
             fprintf(fp,"  QUALIFIER ORIGIN_1");
           else
             fprintf(fp,"  QUALIFIER ORIGIN_2");
           }
      }
    fprintf(fp,"\n");
    }
  else
    fprintf(fp,"%f  %f  %f\n", data->pbc[0], data->pbc[1], R2D*data->pbc[5]);
  }
else
  {
  fputs("CELL\n",fp);
  fprintf(fp,"  %f   %f   %f   90.00000   90.00000   90.00000\n",CELL_LENGTH, CELL_LENGTH, CELL_LENGTH);
  fputs("SYMMETRY  LABEL  P1\n",fp);
  }
 
/* coords */
fprintf(fp,"\nATOMS\n");
if( data->periodic > 0 )
   fprintf(fp,"NAME       X          Y          Z     CHARGE   TEMP    OCCUP   SCAT\n");
else
   fprintf(fp,"NAME      CARX       CARY       CARZ   CHARGE   TEMP    OCCUP   SCAT\n");

/* Cerius2 expects the opposite coordinate ordering for surfaces */
/* is it a surface with negative z? */
if (data->periodic == 2)
  {
  if (g_slist_length(data->cores))
    {
    core = g_slist_nth_data(data->cores, 0);
    if (core->x[2] < 0.0)
      invert=TRUE;
    }
  }

/* rot to make z positive */
matrix_rotation(&tmat[0], PI, PITCH);

for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

  if (core->status & DELETED)
    continue;
  if (core->primary)
    {
/* retrieve */
    ARR3SET(vec, core->x);
    if (invert)
      vecmat(tmat, vec);

    if( !data->periodic )
      { /* Shift non-periodic structure to centre of cell */
      vec[0] += CELL_LENGTH / 2.0;
      vec[1] += CELL_LENGTH / 2.0;
      vec[2] += CELL_LENGTH / 2.0;
      }

    get_elem_data(core->atom_code, &elem_data, data);
    symbol = g_strdup(elem_data.symbol);
    charge = (gint)abs(elem_data.charge);

    fprintf(fp,"%4s %10.5f %10.5f %10.5f %7.4f  0.0000  1.0000   %-2s%d%c\n",
      core->atom_label, vec[0], vec[1], vec[2], atom_charge(core), symbol,
        charge, (elem_data.charge>=0?'+':'-'));
    }
  }
fprintf(fp,"EOF\n");

/* done */
fclose(fp);
return(0);
}

/****************/
/* read routine */
/****************/
#define DEBUG_READ_XTL 0
gint read_xtl(gchar *filename, struct model_pak *data)
{
gint i, offset, num_tokens, num_columns, n=0;
guint shift;
gchar **buff;
gchar *symbol;
struct core_pak *core;
FILE *fp;
gint column[9];

/* initialize columns */
for(i=9;i--;)
  column[i] = -1;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 1);

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* defaults */
data->fractional = TRUE;
data->periodic = 3;

for (;;)
  {
  buff = get_tokenized_line(fp, &num_tokens);

  if (!num_tokens)
    break;

  if (g_ascii_strncasecmp("tit", *buff, 3) == 0) /* TITLE keyword */
    {
    data->title = g_strdup(*(buff+1));
    for(i=2; i<num_tokens; i++)
      data->title = g_strjoin(" ", data->title, *(buff+i), NULL);
    }
/* periodicity search */
  if (g_ascii_strncasecmp("dim", *buff, 3) == 0) /* DIMENSION keyword */
    {
    if (num_tokens > 1)
      data->periodic = (gint) str_to_float(*(buff+1));
    }

/* space group search */
  if (g_ascii_strncasecmp("sym", *buff, 3) == 0) /* SYMMETRY keyword */
    {
    for (i=1 ; i<num_tokens-1 ; i++)
      {
/* seek space group number */
      if (g_ascii_strncasecmp("num", *(buff+i), 3) == 0) /* NUMBER secondary keyword */
        data->sginfo.spacenum = (gint) str_to_float(*(buff+i+1));
/* seek space group label */
      if (g_ascii_strncasecmp("lab", *(buff+i), 3) == 0) /* LABEL secondary keyword */
        data->sginfo.spacename = g_strdup(*(buff+i+1));

/* If space group name not given, use number */
/*      if( strcmp(data->sginfo.spacename,"P 1") == 0 && data->sginfo.spacenum > 1 )
        {
        g_free(data->sginfo.spacename);
        data->sginfo.spacename = g_strdup_printf("%d",data->sginfo.spacenum);
        } */

/* Cell choice options added by C. Fisher 2004 */
/* TODO Read SYM MAT data */
      if (g_ascii_strncasecmp("qua", *(buff+i), 3) == 0) /* QUALIFIER secondary keyword */
        {
        if(g_ascii_strncasecmp("origin_", *(buff+i+1), 7) == 0)
          {
          data->sginfo.cellchoice = (gint) str_to_float(*(buff+i+1)+7);
          sprintf(data->sginfo.spacename, "%s:%d", data->sginfo.spacename, data->sginfo.cellchoice);
          }

        if(g_ascii_strcasecmp("hexagonal", *(buff+i+1)) == 0)
          {
          data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":h", NULL);
          data->sginfo.cellchoice = 1;
          }

        if(g_ascii_strcasecmp("rhombohedral", *(buff+i+1)) == 0)
          {
          data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":r", NULL);
          data->sginfo.cellchoice = 2;
          }

        if(g_ascii_strcasecmp("b_unique", *(buff+i+1)) == 0)
          {
          if( data->sginfo.spacenum == 3 || data->sginfo.spacenum == 4 ||
              data->sginfo.spacenum == 6 || data->sginfo.spacenum == 10 ||
              data->sginfo.spacenum == 11 )
                {
                data->sginfo.cellchoice = 1;
                }
          /* Options for SG 5,7,8,9,12,13,14,15 */
          else if( g_ascii_strncasecmp("c2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("c121",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pc",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p1c1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("cm",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("c1m1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("cc",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("c1c1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c12/m1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/c",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p12/c1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/c",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p121/c1",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("c2/c",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c12/c1",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 1;
                }
          else if( g_ascii_strncasecmp("a2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a121",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pn",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p1n1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("am",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a1m1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("an",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a1n1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a12/m1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p12/n1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/n",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p121/n1",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("a2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a12/n1",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 2;
                }
          else if( g_ascii_strncasecmp("i2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i121",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pa",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p1a1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("im",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i1m1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("ia",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i1a1",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i12/m1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p12/a1",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/a",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p121/a1",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("i2/a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i12/a1",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 3;
                }
              else
                data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":b", NULL);
          }
        if(g_ascii_strcasecmp("c_unique", *(buff+i+1)) == 0)
          {
          if( data->sginfo.spacenum == 3 || data->sginfo.spacenum == 4 ||
              data->sginfo.spacenum == 6 || data->sginfo.spacenum == 10 ||
              data->sginfo.spacenum == 11 )
                {
                data->sginfo.cellchoice = 2;
                }
          /* Options for SG 5,7,8,9,12,13,14,15 */
          else if( g_ascii_strncasecmp("a2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a112",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pa",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p11a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("am",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a11m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("aa",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("a11a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a112/m",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p112/a",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/a",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p1121/a",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("a2/a",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("a112/a",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 4;
                }
          else if( g_ascii_strncasecmp("b2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("b112",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pn",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p11n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("bm",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("b11m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("bn",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("b11n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b112/m",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p112/n",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/n",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p1121/n",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("b2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b112/n",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 5;
                }
          else if( g_ascii_strncasecmp("i2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i112",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pb",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("p11b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("im",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i11m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("ib",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i11b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i112/m",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p112/b",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/b",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p1121/b",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("i2/b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i112/b",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 6;
                }
              else
                data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":c", NULL);
          }
        if(g_ascii_strcasecmp("a_unique", *(buff+i+1)) == 0)
          {
          if( data->sginfo.spacenum == 3 || data->sginfo.spacenum == 4 ||
              data->sginfo.spacenum == 6 || data->sginfo.spacenum == 10 ||
              data->sginfo.spacenum == 11 )
                {
                data->sginfo.cellchoice = 3;
                }
          /* Options for SG 5,7,8,9,12,13,14.15 */
          else if( g_ascii_strncasecmp("b2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("b211",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pb",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("pb11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("bm",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("bm11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("bb",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("bb11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b2/m11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p2/b11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/c",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p21/c11",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("b2/b",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("b2/b11",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 7;
                }
          else if( g_ascii_strncasecmp("c2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("c211",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pn",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("pn11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("cm",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("cm11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("cn",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("cn11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c2/m11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p2/n11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/n",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p21/n11",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("c2/n",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("c2/n11",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 8;
                }
          else if( g_ascii_strncasecmp("i2",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("i211",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("pc",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("pc11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("im",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("im11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("ic",data->sginfo.spacename,2) == 0 ||
              g_ascii_strncasecmp("ic11",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i2/m",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i2/m11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p2/c",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("p2/c11",data->sginfo.spacename,6) == 0 ||
              g_ascii_strncasecmp("p21/c",data->sginfo.spacename,5) == 0 ||
              g_ascii_strncasecmp("p21/c11",data->sginfo.spacename,7) == 0 ||
              g_ascii_strncasecmp("i2/c",data->sginfo.spacename,4) == 0 ||
              g_ascii_strncasecmp("i2/c11",data->sginfo.spacename,6) == 0 )
                {
                data->sginfo.cellchoice = 9;
                }
              else
                data->sginfo.spacename = g_strconcat(data->sginfo.spacename, ":a", NULL);
          }
        }
      }
    }

  if( data->sginfo.spacenum == 50 || data->sginfo.spacenum == 59 ||
      data->sginfo.spacenum == 68 )
    {
    /* Set cellchoice to baseline 1 or 2 */
    if( data->sginfo.cellchoice % 2 == 1 )
       data->sginfo.cellchoice = 1;
    else
       data->sginfo.cellchoice = 2;
    if( g_ascii_strncasecmp(data->sginfo.spacename,"P N C B",7) == 0 ||
        g_ascii_strncasecmp(data->sginfo.spacename,"P N M M",7) == 0 ||
        g_ascii_strncasecmp(data->sginfo.spacename,"C C C B",7) == 0 )
       data->sginfo.cellchoice += 2;
    if( g_ascii_strncasecmp(data->sginfo.spacename,"P C N A",7) == 0 ||
        g_ascii_strncasecmp(data->sginfo.spacename,"P M N M",7) == 0 ||
        g_ascii_strncasecmp(data->sginfo.spacename,"A B A A",7) == 0 )
       data->sginfo.cellchoice += 4;
    if( g_ascii_strncasecmp(data->sginfo.spacename,"A C A A",7) == 0)
       data->sginfo.cellchoice += 6;
    if( g_ascii_strncasecmp(data->sginfo.spacename,"B B C B",7) == 0)
       data->sginfo.cellchoice += 8;
    if( g_ascii_strncasecmp(data->sginfo.spacename,"B B A B",7) == 0)
       data->sginfo.cellchoice += 10;
    }

/* pbc search */
  if (g_ascii_strncasecmp("cel", *buff, 3) == 0) /* CELL keyword */
    {
/* use this line if it has more tokens, otherwise get the next */
    offset = 1;
    if (num_tokens == 1)
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);
      if (!buff)
        return(2);
      offset = 0;
      }

    if (num_tokens > 1)
      {
      data->pbc[0] = str_to_float(*(buff+offset));
      data->pbc[1] = str_to_float(*(buff+offset+1));
      }
    else
      printf("Insufficient CELL tokens.\n");

    switch(data->periodic)
      {
      case 2:
        if (num_tokens > 2)
          data->pbc[5] = D2R * str_to_float(*(buff+offset+2));
        else
          printf("Insufficient CELL tokens.\n");
        break;
      case 3:
        if (num_tokens > 5)
          {
          data->pbc[2] = str_to_float(*(buff+offset+2));
          data->pbc[3] = D2R * str_to_float(*(buff+offset+3));
          data->pbc[4] = D2R * str_to_float(*(buff+offset+4));
          data->pbc[5] = D2R * str_to_float(*(buff+offset+5));
          }
        else
          printf("Insufficient CELL tokens.\n");
        break;
      }
    }

/* coordinate search */
  if (g_ascii_strncasecmp("ato", *buff, 3) == 0) /* ATOMS keyword */
    {
    g_strfreev(buff);
    buff = get_tokenized_line(fp, &num_tokens);
    /* Determine order of columns */
    for(i=0; i<num_tokens; i++)
      {
      /* NAME secondary keyword */
      if (g_ascii_strncasecmp("nam", *(buff+i), 3) == 0)
         column[0] = i;
      /* X secondary keyword */
      if (g_ascii_strncasecmp("x", *(buff+i), 1) == 0)
         column[1] = i;
      /* Y secondary keyword */
      if (g_ascii_strncasecmp("y", *(buff+i), 1) == 0)
         column[2] = i;
      /* Z secondary keyword */
      if (g_ascii_strncasecmp("z", *(buff+i), 1) == 0)
         column[3] = i;
      /* SCATTER secondary keyword */
      if (g_ascii_strncasecmp("sca", *(buff+i), 3) == 0)
         column[4] = i;
      /* CHARGE secondary keyword */
      if (g_ascii_strncasecmp("cha", *(buff+i), 3) == 0)
         column[5] = i;
      /* CARX secondary keyword */
      if (g_ascii_strncasecmp("carx", *(buff+i), 4) == 0)
         {
         column[6] = i;
         data->fractional = FALSE;
         }
      /* CARY secondary keyword */
      if (g_ascii_strncasecmp("cary", *(buff+i), 4) == 0)
         {
         column[7] = i;
         data->fractional = FALSE;
         }
      /* CARZ secondary keyword */
      if (g_ascii_strncasecmp("carz", *(buff+i), 4) == 0)
         {
         column[8] = i;
         data->fractional = FALSE;
         }
      }
    for (;;)
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_columns);
      if (!buff)
        break;
      if (g_ascii_strncasecmp("eof", *buff, 3) == 0)
        break;

/* enough tokens */
/* Use scattering element in lookup, as label not always element symbol */
/* otherwise check if first item is a valid atom type */
/* C. Fisher 2004 */
      if (column[4] != -1)
        {
        symbol = g_strdup(*(buff+column[4]));
        if( num_columns > num_tokens)
          symbol = g_strconcat(symbol, *(buff+column[4]+1), NULL);
        }
      else
        if(column[0] != -1)
          symbol = g_strdup(*(buff+column[0]));
        else
          break;

      if( num_columns > num_tokens) /* SCATTER factor can exist as two items */
        shift=1;      /* Shift column by one to allow for extra scatter term */
      else
        shift=0;

/* add new atom if first item is a valid atom type */
      if (elem_symbol_test(symbol))
        {
        core = new_core(symbol, data);
        data->cores = g_slist_prepend(data->cores, core);

        if( data->fractional )
          {
          if( column[1] > column[4] && shift==1 )
            core->x[0] = str_to_float(*(buff+column[1]+shift));
          else
            core->x[0] = str_to_float(*(buff+column[1]));
          if( column[2] > column[4] && shift==1 )
            core->x[1] = str_to_float(*(buff+column[2]+shift));
          else
            core->x[1] = str_to_float(*(buff+column[2]));
          if( column[3] > column[4] && shift==1 )
            core->x[2] = str_to_float(*(buff+column[3]+shift));
          else
            core->x[2] = str_to_float(*(buff+column[3]));
          }
        else
          {
          if( column[6] > column[4] && shift==1 )
            core->x[0] = str_to_float(*(buff+column[6]+shift));
          else
            core->x[0] = str_to_float(*(buff+column[6]));
          if( column[7] > column[4] && shift==1 )
            core->x[1] = str_to_float(*(buff+column[7]+shift));
          else
            core->x[1] = str_to_float(*(buff+column[7]));
          if( column[8] > column[4] && shift==1 )
            core->x[2] = str_to_float(*(buff+column[8]+shift));
          else
            core->x[2] = str_to_float(*(buff+column[8]));
          }
        if( column[5] != -1)
          {
          core->lookup_charge = FALSE;
          if( column[5] > column[4] && shift==1 )
            core->charge = str_to_float(*(buff+column[5]+shift));
          else
            core->charge = str_to_float(*(buff+column[5]));
          }
        n++;
        }
      }
    }
  }

/* got everything */
data->num_asym = data->num_atoms = n;

#if DEBUG_READ_XTL
printf("SGname %s SGnum %d Option %d\n", data->sginfo.spacename,
  data->sginfo.spacenum, data->sginfo.cellchoice);
printf("Found %d atoms.\n", n);
#endif

/* model setup */
strcpy(data->filename, filename);
g_free(data->basename);
data->basename = parse_strip(filename);
model_prep(data);

return(0);
}
