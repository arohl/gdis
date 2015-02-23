/*
Copyright (C) 2000 by Sean David Fleming

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
#include <unistd.h>

#include "gdis.h"
#include "interface.h"

/* top level data structure */
extern struct sysenv_pak sysenv;

gint error_table_display = TRUE;
GHashTable *error_table=NULL;

#define MAX_ERROR_COUNT 1000

/********************************/
/* error table printing display */
/********************************/
void error_table_disable(void)
{
error_table_display = FALSE;
}

/********************************/
/* error table printing display */
/********************************/
void error_table_enable(void)
{
error_table_display = TRUE;
}

/****************************************/
/* remove all entries in the hash table */
/****************************************/
void error_table_clear(void)
{
if (error_table)
  g_hash_table_destroy(error_table);
error_table=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
}

/************************/
/* print a single entry */
/************************/
void error_table_print(gpointer key, gpointer value, gpointer data)
{
gchar *msg, *text = key;
gint *count = value;

/* TODO - strstr - if error -> gui_show_text ERROR/WARNING/etc */

if (*count > MAX_ERROR_COUNT)
  msg = g_strdup_printf("[count >%d] %s", MAX_ERROR_COUNT, text);
else
  msg = g_strdup_printf("[count %d] %s", *count, text);

gui_text_show(ERROR,msg);

g_free(msg);
}

/******************************************/
/* display all entries in the error table */
/******************************************/
void error_table_print_all(void)
{
/* allow display to be disabled (eg animations) */
if (error_table_display)
  g_hash_table_foreach(error_table, error_table_print, NULL);
}

/**************************************/
/* add a new entry to the error table */
/**************************************/
void error_table_entry(const gchar *text)
{
gint *count;

count = g_hash_table_lookup(error_table, text);
if (!count)
  {
  count = g_malloc0(sizeof(gint));
  g_hash_table_insert(error_table, g_strdup(text), count);
  }
(*count)++;
}

