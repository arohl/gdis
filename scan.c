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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "gdis.h"
#include "file.h"
#include "parse.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***************************/
/* file scanning structure */
/***************************/
struct scan_pak 
{
FILE *fp;

/* TODO - IF decide to do background reads - will need a mutex lock as well */
gboolean background;
gboolean eof;

guint bytes_read;
guint bytes_total;

/* stored lines */
gint buffer_line;
gint buffer_size;
gint buffer_max;
gchar **buffer;
GSList *offset_list;
};

/***********************/
/* setup a new scanner */
/***********************/
gpointer scan_new(gchar *filename)
{
gint i;
struct scan_pak *scan;
struct stat buff;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(NULL);

if (stat(filename, &buff))
  return(NULL);

scan = g_malloc(sizeof(struct scan_pak));
scan->fp = fp;
scan->background = FALSE;
scan->bytes_read = 0;
scan->bytes_total = buff.st_size;
scan->eof = FALSE;
scan->buffer_line = -1;
scan->buffer_size = 0;
scan->buffer_max = 10;
scan->buffer = g_malloc(scan->buffer_max * sizeof(gchar *));
for (i=scan->buffer_max ; i-- ; )
  scan->buffer[i] = NULL;
scan->offset_list = NULL;

return(scan);
}

/*********************/
/* cleanup a scanner */
/*********************/
void scan_free(gpointer ptr_scan)
{
gint i;
struct scan_pak *scan = ptr_scan;

/* checks */
g_assert(scan != NULL);

/* cleanup */
fclose(scan->fp);
for (i=scan->buffer_size ; i-- ; )
  g_free(scan->buffer[i]);
g_free(scan->buffer);

free_slist(scan->offset_list);

g_free(scan);
}

/***************************/
/* EOF extraction primitve */
/***************************/
gboolean scan_complete(gpointer ptr_scan)
{
struct scan_pak *scan = ptr_scan;

/* checks */
if (!scan)
  return(TRUE);

/* completed if we've got an EOF AND we've scanned to the end of the buffer */
if (scan->eof && scan->buffer_line == scan->buffer_size-1)
  return(TRUE);
return(FALSE);
}

/***************************************/
/* retrieve a pointer to the next line */
/***************************************/
/* TODO - return const type as the return should not be free'd by the caller */
#define DEBUG_SCAN_GET_LINE 0
gchar *scan_get_line(gpointer ptr_scan)
{
gint i;
gchar *line;
fpos_t *offset;
struct scan_pak *scan = ptr_scan;

/* checks */
g_assert(scan != NULL);

/* TODO - progress bar (popup or part of main gdis win) for large files */
/*
printf("read %d/%d\n", scan->bytes_read, scan->bytes_total);
*/

if (scan->buffer_line >= scan->buffer_size)
  {
printf("premature EOF?\n");
g_assert_not_reached();
  }

g_assert(scan->buffer_line < scan->buffer_size);

/* read new line, or can we do a buffered read? */
if (scan->buffer_line < scan->buffer_size-1)
  {
  g_assert(scan->buffer_line >= 0);

/* get next line in the buffer */
  scan->buffer_line++;
  line = scan->buffer[scan->buffer_line];
  }
else
  {
/* read new line */
  if (scan->buffer_size == scan->buffer_max)
    {
/* max size reached; discard oldest line and shuffle the buffer */
    g_free(scan->buffer[0]);
    for (i=1 ; i<scan->buffer_max ; i++)
      scan->buffer[i-1] = scan->buffer[i];

/* remove first element (oldest item read in) */
/*
    offset = (scan->offset_list)->data;
*/

    offset = g_slist_nth_data(scan->offset_list, 0);

    scan->offset_list = g_slist_remove(scan->offset_list, offset);
    g_free(offset);
    }
  else
    {
/* keep filling buffer until max size reached */
    scan->buffer_line++;
    scan->buffer_size++;
    }

/* store file pointer offset for the current line */
  offset = g_malloc(sizeof(fpos_t));
  if (!fgetpos(scan->fp, offset))
    scan->offset_list = g_slist_append(scan->offset_list, offset);
  else
    {
    g_free(offset);
    printf("WARNING: failed to save file pointer offset; animations won't work!\n");
    }

/* read a new line into the buffer */
  line = file_read_line(scan->fp);
  if (!line)
    scan->eof = TRUE;
  else
    scan->bytes_read += strlen(line);

  g_assert(scan->buffer_line >= 0);

  scan->buffer[scan->buffer_line] = line;
  }

#if DEBUG_SCAN_GET_LINE 
printf("> %s", line);
#endif
return(line);
}

/*************************************/
/* get next (non white space) tokens */
/*************************************/
gchar **scan_get_tokens(gpointer ptr_scan, gint *num_tokens)
{
gchar *line, **buff;
struct scan_pak *scan = ptr_scan;

g_assert(scan != NULL);

while (!scan_complete(scan))
  {
  line = scan_get_line(scan);
  buff = tokenize(line, num_tokens);
  if (buff)
    return(buff);
  }
*num_tokens = 0;
return(NULL);
}

/*******************************/
/* get pointer to current line */
/*******************************/
gchar *scan_cur_line(gpointer ptr_scan)
{
struct scan_pak *scan = ptr_scan;

g_assert(scan != NULL);

return(scan->buffer[scan->buffer_line]);
}

/******************************************/
/* flag the current line as a frame start */
/******************************************/
void scan_frame_new(gpointer ptr_scan, struct model_pak *model)
{
gpointer offset, data;
struct scan_pak *scan = ptr_scan;

g_assert(scan != NULL);
g_assert(model != NULL);

data = g_slist_nth_data(scan->offset_list, scan->buffer_line);

offset = g_memdup(data, sizeof(fpos_t));

model->frame_list = g_list_append(model->frame_list, offset); 

/*
printf("Add offset: %p [%d/%d]\n", offset, scan->buffer_line, g_slist_length(scan->offset_list));
*/
}

/*************************/
/* store a file position */
/*************************/
gpointer scan_offset_get(gpointer ptr_scan)
{
struct scan_pak *scan = ptr_scan;

printf("buffer line: %d\n", scan->buffer_line);
printf("buffer length: %d\n", scan->buffer_size);
printf("offset length: %d\n", g_slist_length(scan->offset_list));

return(g_slist_nth_data(scan->offset_list, scan->buffer_line));
}

/*************************************/
/* put back a line (into the buffer) */
/*************************************/
gboolean scan_put_line(gpointer ptr_scan)
{
struct scan_pak *scan = ptr_scan;

g_assert(scan != NULL);

if (scan->buffer_line > 0)
  scan->buffer_line--;
else
  return(TRUE);

return(FALSE);
}

