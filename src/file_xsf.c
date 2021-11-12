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
#include <string.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "matrix.h"
#include "interface.h"
#include "parse.h"
#include "space.h"
#include "edit.h"

/* main structures */
extern struct sysenv_pak sysenv;

/****************************/
/* write a coordinate block */
/****************************/
gint write_xsf_frame(FILE *fp, struct model_pak *model)
{
gint i=0;
gdouble vec[3];
GSList *list;
struct core_pak *core;

/* header */
if (model->periodic == 3)
  {
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = (struct core_pak *) list->data;

    if (core->status & (DELETED | HIDDEN))
      continue;
    i++;
    }
  fprintf(fp,"%4d 1\n", i);
  }

/* coords */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;

  if (core->status & (DELETED | HIDDEN))
    continue;
  if (core->primary)
    {
    ARR3SET(vec, core->x);
    vecmat(model->latmat, vec);

    fprintf(fp,"%4d %9.5f %9.5f %9.5f\n",
            core->atom_code, vec[0], vec[1], vec[2]);
    }
  }
return(0);
}

/**********************/
/* save in XSF format */
/**********************/
#define DEBUG_WRITE_XSF 0
gint write_xsf(gchar *filename, struct model_pak *data)
{
gint n, current, status;
FILE *fp;
GString *err_text;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

if( g_slist_length(data->cores) > 9999)
  {
  gui_text_show(ERROR, "Too many atoms for XSF format!\n");
  return(1);
  }

if( data->periodic == 1)
  {
  gui_text_show(ERROR, "Model must be 0, 2 or 3 dimensional.\n");
  return(1);
  }

/* open the file */
fp = fopen(filename, "w");
if (!fp)
  {
  gui_text_show(ERROR, "Bad filename!\n");
  return(1);
  }

if (data->frame_list)
  {
  data->afp = fopen(data->filename, "r");

  if (!data->afp)
    {
    gui_text_show(ERROR, "Failed to open animation stream.\n");
    return(0);
    }

  current = data->cur_frame; /* Remember current position */

  fprintf(fp, "ANIMSTEPS %d\n", data->num_frames);

  for (n = data->cur_frame ; n < data->num_frames ; n++)
    {
    status = read_raw_frame(data->afp, n, data);
#if DEBUG_WRITE_XSF
printf("Frame no: %d\n",n);
#endif
    if (status)
      {
      err_text = g_string_new("");
      g_string_printf(err_text, "Error reading frame: %d\n", n);
      gui_text_show(ERROR, err_text->str);
      g_string_free(err_text, TRUE);
      return(2);
      }

    if (data->periodic == 0)
      {
      fprintf(fp, "ATOMS %d\n", n+1);
      }
    else
      {
      if (data->periodic == 2)
         fputs("SLAB\n", fp);
      else
         fputs("CRYSTAL\n", fp);

      fprintf(fp, "PRIMVEC %d\n", n+1);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[0], data->latmat[3], data->latmat[6]);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[1], data->latmat[4], data->latmat[6]);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[2], data->latmat[5], data->latmat[8]);

      if (data->periodic == 3)
        {
        fprintf(fp, "CONVVEC %d\n", n+1);
        fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[0], data->latmat[3], data->latmat[6]);
        fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[1], data->latmat[4], data->latmat[7]);
        fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[2], data->latmat[5], data->latmat[8]);
        }

      fprintf(fp, "PRIMCOORD %d\n", n+1);
      }

    write_xsf_frame(fp, data);
    }

  /* Restore data for current frame */
  read_raw_frame(data->afp, current, data);
  connect_molecules(data);
  fclose(data->afp);
  }
else
  {
  if (data->periodic == 0)
    fputs("ATOMS\n", fp);
  else
    {
    if (data->periodic == 2)
      fputs("SLAB\n", fp);
    else
      fputs("CRYSTAL\n", fp);

    fputs("PRIMVEC\n", fp);
    fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[0], data->latmat[3], data->latmat[6]);
    fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[1], data->latmat[4], data->latmat[7]);
    fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[2], data->latmat[5], data->latmat[8]);

    if (data->periodic == 3)
      {
      fputs("CONVVEC\n", fp);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[0], data->latmat[3], data->latmat[6]);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[1], data->latmat[4], data->latmat[7]);
      fprintf(fp, "%11.6f %11.6f %11.6f\n", data->latmat[2], data->latmat[5], data->latmat[8]);
      }
    fputs("PRIMCOORD\n", fp);
    }
  write_xsf_frame(fp, data);
  }

/* done */
gui_text_show(STANDARD, "File saved in XSF format.\n");
fclose(fp);
return(0);
}

/*********************************************/
/* read a xsf block into the model structure */
/*********************************************/
#define DEBUG_READ_XSF_BLOCK 0
/* NB: assumes fp is one line before frame data */
void read_xsf_block(FILE *fp, struct model_pak *data)
{
gint i, j, k, code, num_tokens;
gint n = 0, num_atoms = 0;
//gint frame = 0;
gint convflag = 0; /* Conventional cell given = 1 */
gchar **buff;
gdouble m[3];
gdouble convmat[9], invmat[9];
struct core_pak *core, *copy;
struct elem_pak elem_data;
GSList *list, *ilist = NULL;

g_assert(fp != NULL);
g_assert(data != NULL);

if (data->cores)
  free_core_list(data);

/* loop while there's data */
for (;;)
  {
  buff = get_tokenized_line(fp, &num_tokens);
  if (!buff)
    break;

  if (num_tokens)
    {
    if (!g_ascii_strcasecmp(*buff, "atoms") && n > 0)
      {
      fseek(fp, ftell(fp)-1, SEEK_SET); /* rewind a line to preserve tag info */
      g_strfreev(buff);
      break;
      }

    if (!g_ascii_strcasecmp(*buff, "primvec"))
      {
      data->construct_pbc = TRUE;
      for (i=0; i<3; i++)
        {
        g_strfreev(buff);
        buff = get_tokenized_line(fp, &num_tokens);
        if (num_tokens == 3)
          {
          data->latmat[i] = str_to_float(*(buff));
          data->latmat[i+3] = str_to_float(*(buff+1));
          data->latmat[i+6] = str_to_float(*(buff+2));
          }
        }
#if DEBUG_READ_XSF_BLOCK
P3MAT("Primitive lattice:", data->latmat);
#endif
      continue;
      }

    if (!g_ascii_strcasecmp(*buff, "convvec"))
      {
      convflag = 1;
      for (i=0; i<3; i++)
        {
        g_strfreev(buff);
        buff = get_tokenized_line(fp, &num_tokens);
        if (num_tokens == 3)
          {
          convmat[i] = str_to_float(*(buff));
          convmat[i+3] = str_to_float(*(buff+1));
          convmat[i+6] = str_to_float(*(buff+2));
          }
        }
#if DEBUG_READ_XSF_BLOCK
P3MAT("Conventional lattice:", convmat);
#endif
      continue;
      }
 
    if (!g_ascii_strcasecmp(*buff, "primcoord"))
      {
      g_strfreev(buff);
      buff = get_tokenized_line(fp, &num_tokens);

      if (num_tokens == 2)
        num_atoms = (gint)str_to_float(*buff);
      else
	{
        g_strfreev(buff);
        return;
	}
      continue;
      }

    if (num_tokens > 3)
      {
      code = elem_symbol_test(*buff);

      if (!code)
        code = (gint)str_to_float(*buff);

      if (code > 0)
        {
        n++;
/* otherwise, add a new core */
        get_elem_data(code, &elem_data, data);
        core = core_new(elem_data.symbol, NULL, data);
        data->cores = g_slist_append(data->cores, core);

        core->x[0] = str_to_float(*(buff+1));
        core->x[1] = str_to_float(*(buff+2));
        core->x[2] = str_to_float(*(buff+3));
        core->lookup_charge = TRUE;
#if DEBUG_READ_XSF_BLOCK
printf("Atom %d %s %f %f %f\n", n, core->atom_label,
           core->x[0], core->x[1], core->x[2]);
#endif

/* end this frame if we've got the number of atoms expected */
        if (n == num_atoms)
          {
          g_strfreev(buff);
          break;
          }
        }
      }
    }
  g_strfreev(buff);
  }

#if DEBUG_READ_XSF_BLOCK
printf("Frame: %d\n", data->cur_frame);
#endif

if (n != g_slist_length(data->cores) || 
    (data->periodic > 0 && num_atoms != n) ||
    n == 0)
  gui_text_show(ERROR, "Wrong number of atoms in xsf file.\n");

/* convert input cartesian coords to fractional */
if( data->periodic)
  {
  matrix_lattice_init(data);
  coords_make_fractional(data);
  data->fractional = TRUE;
  }

/* Convert to conventional cell if matrix given */
if (convflag)
  {
  memcpy(invmat, data->latmat, 9*sizeof(gdouble));
  matrix_invert(invmat);
  memcpy(data->ilatmat, invmat, 9*sizeof(gdouble));
  matmat(convmat, invmat);
  for (i=0; i<3; i++)
    m[i] = (gint) ceil(sqrt(invmat[i]+invmat[i+3]+invmat[i+6]))+1;
  /* create image cores */
  for (i = 0; i < m[0]; i++)
    for (j = 0; j < m[1]; j++)
      for (k = 0; k < m[2]; k++)
        {
	if (i == 0 && j == 0 && k == 0)
	  continue;
	else
          for (list=data->cores ; list ; list=g_slist_next(list))
            {
            core = list->data;
            copy = dup_core(core);
            copy->x[0] += i;
            copy->x[1] += j;
            copy->x[2] += k;
            ilist = g_slist_prepend(ilist, copy);
            }
        }

  /* order preserving core/shell additions */
  data->cores = g_slist_concat(data->cores, g_slist_reverse(ilist));

  /* convert to cartesian */
  coords_make_cartesian(data);

  data->periodic = 3;

  /* set latmat to conventional cell values */
  ARR3SET(data->latmat, convmat);
  ARR3SET(&(data->latmat[3]), &convmat[3]);
  ARR3SET(&(data->latmat[6]), &convmat[6]);
  matrix_lattice_init(data);
  coords_make_fractional(data);
  data->fractional = TRUE;
  }

model_prep(data);
}

/******************/
/* xsf frame read */
/******************/
gint read_xsf_frame(FILE *fp, struct model_pak *data)
{
/* frame overwrite */
read_xsf_block(fp, data);
return(0);
}
/****************/
/* read routine */
/****************/
#define DEBUG_READ_XSF 0
gint read_xsf(gchar *filename, struct model_pak *data)
{
gint frame = 0, end_flag = 0;
gint i, total_frames = 1;
gint prim_flag = 0, num_tokens;
gchar **buff;
long int fp_pos, fp_pos_prev = 0;
fpos_t *offset;
FILE *fp;

/* checks */
g_return_val_if_fail(data != NULL, 1);
g_return_val_if_fail(filename != NULL, 1);

if (strlen(filename) < FILELEN)
  strcpy(data->filename, filename);
else
  {
  gui_text_show(ERROR, "File name is too long.\n");
  return(-1);
  }

fp = fopen(filename, "rt");
if (fp == NULL)
  return(2);

/* inits */
data->fractional = FALSE;
data->periodic = 0;

for (;;)
  {
  fp_pos = ftell(fp); /* file position before reading line */
  buff = get_tokenized_line(fp, &num_tokens);
  if (buff == NULL)
    break;

  /* first line should be type */
  if (!g_ascii_strcasecmp(*buff, "animsteps"))
    {
    total_frames = str_to_float(*(buff+1));
    continue;
    }
  else if (!g_ascii_strcasecmp(*buff, "slab"))
    {
    data->periodic = 2;
    continue;
    }
  else if (!g_ascii_strcasecmp(*buff, "crystal"))
    {
    data->periodic = 3;
    continue;
    }
  else if (!g_ascii_strcasecmp(*buff, "atoms"))
    {
    fp_pos_prev = fp_pos; /* file position before reading line */
    end_flag++;
    }
  else if (!g_ascii_strcasecmp(*buff, "primvec"))
    {
    if (num_tokens == 2)
      {
      fp_pos_prev = fp_pos; /* file position before reading line */
      prim_flag++;
      end_flag++;
      }
    else
      {
      data->construct_pbc = TRUE;
      for (i=0; i<3; i++)
        {
        g_strfreev(buff);
        buff = get_tokenized_line(fp, &num_tokens);
        if (num_tokens == 3)
          {
          data->latmat[3*i] = str_to_float(*(buff));
          data->latmat[3*i+1] = str_to_float(*(buff+1));
          data->latmat[3*i+2] = str_to_float(*(buff+2));
          }
        }
#if DEBUG_READ_XSF
P3MAT("Fixed lattice:", data->latmat);
#endif
      continue;
      }
    }
  else if (!g_ascii_strcasecmp(*buff, "primcoord"))
    {
    if (prim_flag == 0)
      {
      fp_pos_prev = fp_pos; /* file position before reading line */
      end_flag++;
      }
    }
  else
    end_flag = 0;

/* coords search */
  if ( end_flag == 1)
    {
    fp_pos = ftell(fp);
    fseek(fp, fp_pos_prev, SEEK_SET); /* rewind a line to preserve tag info */
    add_frame_offset(fp, data);
    fseek(fp, fp_pos, SEEK_SET); /* return to current position */

/* increment counter */
    frame++;
    }

  g_strfreev(buff);
  }

offset = g_list_nth_data(data->frame_list, 0);

#if DEBUG_READ_XSF
printf("Periodicity: %d\n", data->periodic);
printf("Total frames: %d\n", g_list_length(data->frame_list));
#endif

if (!offset)
  return(1);

if (fsetpos(fp, offset))
  {
  printf("Error positioning file pointer.\n");
  return(2);
  }
else
  {
/* use supplied routine (if available) */
  read_xsf_frame(fp, data);
  }

/* got everything */
if (frame != total_frames)
  {
  gui_text_show(ERROR, "Incorrect number of frames.\n");
  return(2);
  }

data->num_asym = g_slist_length(data->cores);
data->num_frames = frame;

/* get rid of frame list if only one frame */
if (data->num_frames == 1)
  {
  free_list(data->frame_list);
  data->frame_list = NULL;
  }

#if DEBUG_READ_XSF
printf("Found %d atoms.\n", g_slist_length(data->cores));
#endif

/* model setup */
g_free(data->basename);
data->basename = parse_strip(filename);

model_prep(data);

return(0);
}
