/*
Copyright (C) 2003 by Andrew Lloyd Rohl

andrew@ivec.org

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

#include "gdis.h"
#include "gamess.h"
#include "coords.h"
#include "model.h"
#include "edit.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "interface.h"


#define GMS_NGAUSS_TXT "ngauss="
#define GMS_NUM_P_TXT "npfunc="
#define GMS_NUM_D_TXT "ndfunc="
#define GMS_NUM_F_TXT "nffunc="
#define GMS_DIFFSP_TXT "diffsp=.true."
#define GMS_DIFFS_TXT "diffs=.true."
/* TODO: POLAR and split */

#define GMS_MPLEVEL_TXT "mplevl="
#define GMS_CITYP_TXT "cityp="
#define GMS_CCTYP_TXT "cctyp="
#define GMS_MAXIT_TXT "maxit="
#define GMS_TOTAL_Q_TXT "icharg="
#define GMS_MULT_TXT "mult="
#define GMS_WIDE_OUTPUT_TXT "nprint="
#define GMS_COORD_TXT "coord=cart"

#define GMS_TIMLIM_TXT "timlim="
#define GMS_MWORDS_TXT "mwords="

#define GMS_NSTEP_TXT "nstep="

#define GMS_UNITS_TXT "units="
struct GMS_keyword_pak units[] = {
{"angs", GMS_ANGS},
{"bohr", GMS_BOHR},
{NULL, 0}
};

#define GMS_EXETYPE_TXT "exetyp="
struct GMS_keyword_pak exe_types[] = {
{"run", GMS_RUN},
{"check", GMS_CHECK},
{"debug", GMS_DEBUG},
{NULL, 0}
};

#define GMS_RUNTYPE_TXT "runtyp="
struct GMS_keyword_pak run_types[] = {
{"energy",   GMS_ENERGY},
{"gradient", GMS_GRADIENT},
{"hessian",  GMS_HESSIAN},
{"prop",     GMS_PROP},
{"morokuma", GMS_MOROKUMA},
{"transitn", GMS_TRANSITN},
{"ffield",   GMS_FFIELD},
{"tdhf",     GMS_TDHF},
{"makefp",   GMS_MAKEFP},
{"optimize", GMS_OPTIMIZE},
{"trudge",   GMS_TRUDGE},
{"sadpoint", GMS_SADPOINT},
{"irc",      GMS_IRC},
{"vscf",     GMS_VSCF},
{"drc",      GMS_DRC},
{"globop",   GMS_GLOBOP},
{"gradextr", GMS_GRADEXTR},
{"surface",  GMS_SURFACE},
{NULL, 0} 
};

#define GMS_SCFTYPE_TXT "scftyp="
struct GMS_keyword_pak scf_types[] = {
{"rhf",  GMS_RHF},
{"uhf",  GMS_UHF}, 
{"rohf", GMS_ROHF}, 
{NULL, 0} 
};

#define GMS_METHOD_TXT "method="
struct GMS_keyword_pak method_types[] = {
{"qa",       GMS_QA}, 
{"nr",       GMS_NR},
{"rfo",      GMS_RFO},
{"schlegel", GMS_SCHLEGEL},
{NULL, 0} 
};

#define GMS_BASIS_TXT "gbasis="
struct GMS_keyword_pak basis_types[] = {
{"user defined", GMS_USER},
{"mndo", GMS_MNDO},
{"am1",  GMS_AM1},
{"pm3",  GMS_PM3},
{"mini", GMS_MINI},
{"midi", GMS_MIDI},
{"sto",  GMS_STO},
{"n21",  GMS_N21},
{"n31",  GMS_N31},
{"n311", GMS_N311},
{"dzv",  GMS_DZV},
{"dh",   GMS_DH},
{"tzv",  GMS_TZV},
{"mc",   GMS_MC},
{NULL, 0} 
};

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* structures for gamess */
struct basis_pak basis_sets[] = {
{"User Defined", GMS_USER, 0},
{"MNDO", GMS_MNDO, 0},
{"AM1", GMS_AM1, 0},
{"PM3", GMS_PM3, 0},
{"MINI", GMS_MINI, 0},
{"MIDI", GMS_MIDI, 0},
{"STO-2G", GMS_STO, 2},
{"STO-3G", GMS_STO, 3},
{"STO-4G", GMS_STO, 4},
{"STO-5G", GMS_STO, 5},
{"STO-6G", GMS_STO, 6},
{"3-21G", GMS_N21, 3},
{"6-21G", GMS_N21, 6},
{"4-31G", GMS_N31, 4},
{"5-31G", GMS_N31, 5},
{"6-31G", GMS_N31, 6},
{"6-311G", GMS_N311, 6},
{"DZV", GMS_DZV, 0},
{"DH", GMS_DH, 0},
{"TZV", GMS_TZV, 0},
{"MC", GMS_MC, 0},
{NULL, 0, 0}
};

/***********************************/
/* GAMESS structure initialization */
/***********************************/
void gamess_init(gpointer data)
{
struct gamess_pak *gamess = data;

gamess->units = GMS_ANGS;
gamess->exe_type = GMS_RUN;
gamess->run_type = GMS_ENERGY;
gamess->scf_type = GMS_RHF;
gamess->basis = GMS_MNDO;
gamess->dft = FALSE;
gamess->dft_functional = 0;
gamess->ngauss = 0;
gamess->opt_type = GMS_QA;
gamess->total_charge = 0;
gamess->multiplicity = 1;
gamess->time_limit = 600;
gamess->mwords = 1;
gamess->num_p = 0;
gamess->num_d = 0;
gamess->num_f = 0;
gamess->have_heavy_diffuse = FALSE;
gamess->have_hydrogen_diffuse = FALSE;
gamess->converged = FALSE;
gamess->wide_output = FALSE;
gamess->nstep = 20;
gamess->maxit = 30;
gamess->title = g_strdup("none");
gamess->temp_file = g_strdup("none");
gamess->out_file = g_strdup("none");
gamess->energy = gamess->max_grad = gamess->rms_grad = 0.0;
gamess->have_energy = gamess->have_max_grad = gamess->have_rms_grad = FALSE;
}

/****************/
/* file writing */
/****************/
void write_keyword(FILE *fp, gchar *keyword, gint id, struct GMS_keyword_pak *values)
{
gint i=0;
while (values[i].label)
  {
  if (values[i].id == id)
    fprintf(fp, " %s%s", keyword, values[i].label);
  i++;
  }
}
 
gint write_gms(gchar *filename, struct model_pak *model)
{
gdouble x[3];
GSList *list;
struct core_pak *core;
FILE *fp;

/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);

/* open the file */
fp = fopen(filename,"wb");
if (!fp)
  return(3);

/* TODO - enforce lines < 80 chars? */  
/* print control keywords */
fprintf(fp, " $contrl coord=unique");
write_keyword(fp, GMS_EXETYPE_TXT, model->gamess.exe_type, exe_types);
write_keyword(fp, GMS_SCFTYPE_TXT, model->gamess.scf_type, scf_types);
write_keyword(fp, GMS_RUNTYPE_TXT, model->gamess.run_type, run_types);

/* FIXME - put in to try & keep the lines below 80 chars */
fprintf(fp, "\n ");

write_keyword(fp, GMS_UNITS_TXT, model->gamess.units, units);

fprintf(fp, " %s%d", GMS_MAXIT_TXT, (gint) model->gamess.maxit);
if (model->gamess.total_charge != 0.0)
  fprintf(fp, " %s%d", GMS_TOTAL_Q_TXT, (gint) model->gamess.total_charge);
if (model->gamess.multiplicity > 1)
  fprintf(fp, " %s%d", GMS_MULT_TXT, (gint) model->gamess.multiplicity);
if (model->gamess.wide_output)
  fprintf(fp, " %s6", GMS_WIDE_OUTPUT_TXT);

fprintf(fp, " $end\n");

/* NEW - dft calc */
if (model->gamess.dft)
  {
  fprintf(fp, " $dft");
  switch (model->gamess.dft_functional)
    {
    case SVWN:
      fprintf(fp, " dfttyp=SVWN");
      break; 
    case BLYP:
      fprintf(fp, " dfttyp=BLYP");
      break; 
    case B3LYP:
      fprintf(fp, " dfttyp=B3LYP");
      break; 
    }
  fprintf(fp, " $end\n");
  }

/* TODO: electron correlation stuff */

/* print size and memory */
fprintf(fp, " $system %s%d %s%d $end\x0a", GMS_TIMLIM_TXT, (gint) model->gamess.time_limit, GMS_MWORDS_TXT, (gint) model->gamess.mwords);

/* print optimiser data */
if (model->gamess.run_type >= GMS_OPTIMIZE)
 {
  fprintf(fp, " $statpt %s%d ", GMS_NSTEP_TXT, (gint) model->gamess.nstep);
  write_keyword(fp, GMS_METHOD_TXT, model->gamess.opt_type, method_types);
  fprintf(fp, "$end\n");
  }

/* print basis set if one of the standard ones */
if (model->gamess.basis != GMS_USER)
  {
  fprintf(fp, " $basis ");
  write_keyword(fp, GMS_BASIS_TXT, model->gamess.basis, basis_types);
  if (model->gamess.ngauss)
    fprintf(fp, " %s%d", GMS_NGAUSS_TXT, model->gamess.ngauss);
  if (model->gamess.num_p)
    fprintf(fp, " %s%d", GMS_NUM_P_TXT, (gint) model->gamess.num_p);
  if (model->gamess.num_d)
    fprintf(fp, " %s%d", GMS_NUM_D_TXT, (gint) model->gamess.num_d);
  if (model->gamess.num_f)
    fprintf(fp, " %s%d", GMS_NUM_F_TXT, (gint) model->gamess.num_f);
  if (model->gamess.have_heavy_diffuse)
    fprintf(fp, " %s", GMS_DIFFSP_TXT);
  if (model->gamess.have_hydrogen_diffuse)
    fprintf(fp, " %s", GMS_DIFFS_TXT);
  fprintf(fp, " $end\n");
  }
  
/* print data */
fprintf(fp, " $DATA\n");
/* print data header */
fprintf(fp, "%s\n", model->gamess.title);
/* print symmetry */
fprintf(fp, "c1\n");

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & DELETED)
    continue;

/* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
  if (model->gamess.units == GMS_ANGS)
    fprintf(fp,"%-7s  %d  %14.9f  %14.9f  %14.9f\n",
              elements[core->atom_code].symbol, elements[core->atom_code].number, x[0], x[1], x[2]);
  else
    fprintf(fp,"%-7s  %d  %14.9f  %14.9f  %14.9f\n",
              elements[core->atom_code].symbol, elements[core->atom_code].number, x[0]/BOHR_TO_ANGS, x[1]/BOHR_TO_ANGS, x[2]/BOHR_TO_ANGS);
  }
fprintf(fp, " $END\n");

fclose(fp);
return(0);
}

gchar *get_next_keyword(FILE *fp, gchar *line, gint newline, gint *ret)
{
static int i;
int num_tokens;
gchar *text;
static gchar **buff;
gchar *keyword;

if (newline)
  {
  buff = tokenize(line, &num_tokens);
  i = 0;
  }
while (buff[i] == NULL)
  {
  /* TODO skip comment lines? */

  text = file_read_line(fp);
  if (!text)
    {
    *ret = 2;
    return(NULL);
    }
  else
    {
    /* TODO free array of pointers */
    /* g_strfreev(buff); */
    buff = tokenize(text, &num_tokens);
    i = 0;
    }
  g_free(text);
  }
*ret = 0;
keyword = buff[i];
i++;
return(keyword);
}

GSList *get_gamess_keywords(FILE *fp, gchar *line, gint *ret)
{
gchar *keyword;
GSList *keywords=NULL;

keyword = get_next_keyword(fp, line, TRUE, ret);
while ((g_ascii_strncasecmp(keyword, "$end", 4) != 0) && (*ret == 0))
  {
  keywords = g_slist_append(keywords, keyword);
  keyword = get_next_keyword(fp, line, FALSE, ret);
  }
if (g_ascii_strncasecmp(keyword, "$end", 4) == 0)
  g_free(keyword);
return(keywords);
}

gint read_keyword(gchar *value, struct GMS_keyword_pak *values, gint *id)
{
gint i=0, found=FALSE;
while (values[i].label)
  {
  if (g_ascii_strcasecmp(values[i].label, value) == 0)
    {
    *id = values[i].id;
    found=TRUE;
    }
  i++;
  }
  if (found)
    return(0);
  else
    return(1);
}

gint get_data(FILE *fp, struct model_pak *model, gint have_basis)
{
gchar **buff, *line;
gint num_tokens;
struct core_pak *core;

/* process title line */
line = file_read_line(fp);
if (!line)
  {
  gui_text_show(ERROR, "unexpected end of file reading title\n");
  return(2);
  }
model->gamess.title = g_strstrip(line);

/* process symmetry line */
line = file_read_line(fp);
if (!line)
  {
  gui_text_show(ERROR, "unexpected end of file reading symmetry\n");
  return(2);
  }
if (g_ascii_strncasecmp(line, "c1", 2) != 0)
  {
  /* TODO handle symmetry! */
  gui_text_show(WARNING, "only C1 symmetry understood at present\n");
  }
g_free(line);

/* process coordinates */
line = file_read_line(fp);
if (!line)
  {
  gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
  return(2);
  }

while (g_ascii_strncasecmp(line, " $end", 5) != 0)
  {
  buff = tokenize(line, &num_tokens);
  /* TODO Store GAMESS label and use in .inp files */

  if (num_tokens > 4)
    {
    core = new_core(elements[(int) str_to_float(*(buff+1))].symbol, model);
    if (model->gamess.units == GMS_ANGS)
      {
      core->x[0] = str_to_float(*(buff+2));
      core->x[1] = str_to_float(*(buff+3));
      core->x[2] = str_to_float(*(buff+4));
      }
    else
      {
      core->x[0] = str_to_float(*(buff+2)) * BOHR_TO_ANGS;
      core->x[1] = str_to_float(*(buff+3)) * BOHR_TO_ANGS;
      core->x[2] = str_to_float(*(buff+4)) * BOHR_TO_ANGS;
      }
    model->cores = g_slist_append(model->cores, core);
    g_strfreev(buff);
    if (!have_basis)
      {
      model->gamess.basis = GMS_USER;

/* TODO - read instead of skipping */
      for(;;)
        {
        line = file_read_line(fp);
        if (!line)
          break;

        buff = tokenize(line, &num_tokens);
        g_strfreev(buff);
        if (!num_tokens)
          break;
        }
      }
    }
/* get the next line */
  g_free(line);
  line = file_read_line(fp);
  if (!line)
    {
    gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
    return(2);
    }
  }
g_free(line);
return(0);
}

gint get_basis(gchar *keyword, struct model_pak *model)
{
gint len, basis;

if (g_ascii_strncasecmp(GMS_BASIS_TXT, keyword, len=strlen(GMS_BASIS_TXT)) == 0)
  {
  if (read_keyword(&keyword[len], basis_types, &basis) > 0)
    {
    gui_text_show(ERROR, " unknown basis ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  model->gamess.basis = basis;
  return(0);
  }
else if (g_ascii_strncasecmp(keyword, GMS_NGAUSS_TXT, len=strlen(GMS_NGAUSS_TXT)) == 0)
  model->gamess.ngauss = (gint) (str_to_float(&keyword[len]));
else if (g_ascii_strncasecmp(keyword, GMS_NUM_P_TXT, len=strlen(GMS_NUM_P_TXT)) == 0)
  model->gamess.num_p = str_to_float(&keyword[len]);
else if (g_ascii_strncasecmp(keyword, GMS_NUM_D_TXT, len=strlen(GMS_NUM_D_TXT)) == 0)
  model->gamess.num_d = str_to_float(&keyword[len]);
else if (g_ascii_strncasecmp(keyword, GMS_NUM_F_TXT, len=strlen(GMS_NUM_F_TXT)) == 0)
  model->gamess.num_f = str_to_float(&keyword[len]);
else if (g_ascii_strncasecmp(keyword, GMS_DIFFSP_TXT, strlen(GMS_DIFFSP_TXT)) == 0)//len UNUSED
  model->gamess.have_heavy_diffuse = TRUE;
else if (g_ascii_strncasecmp(keyword, GMS_DIFFS_TXT, strlen(GMS_DIFFS_TXT)) == 0)//len UNUSED
  model->gamess.have_hydrogen_diffuse = TRUE;
else 
  {
  gui_text_show(ERROR, " unknown keyword ");
  gui_text_show(ERROR, keyword);
  gui_text_show(ERROR, " ");
  return(1);
  }
return(0);
}

gint get_control(gchar *keyword, struct model_pak *model)
{
gint len, bad_andrew;

if (g_ascii_strncasecmp(GMS_UNITS_TXT, keyword, len=strlen(GMS_UNITS_TXT)) == 0)
  {
  if (read_keyword(&keyword[len], units, &bad_andrew) > 0)
    {
    gui_text_show(ERROR, " unknown units ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  model->gamess.units = bad_andrew;
  return(0);
  }
if (g_ascii_strncasecmp(GMS_EXETYPE_TXT, keyword, len=strlen(GMS_EXETYPE_TXT)) == 0)
  {
  if (read_keyword(&keyword[len], exe_types, &bad_andrew) > 0)
    {
    gui_text_show(ERROR, " unknown exetyp ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  model->gamess.exe_type = bad_andrew;
  return(0);
  }

if (g_ascii_strncasecmp(GMS_RUNTYPE_TXT, keyword, len=strlen(GMS_RUNTYPE_TXT)) == 0)
  {
  if (read_keyword(&keyword[len], run_types, &bad_andrew) > 0)
    {
    gui_text_show(ERROR, " unknown runtyp ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  model->gamess.run_type = bad_andrew;
  return(0);
  }

if (g_ascii_strncasecmp(GMS_SCFTYPE_TXT, keyword, len=strlen(GMS_SCFTYPE_TXT)) == 0)
  {
  if (read_keyword(&keyword[len], scf_types, &bad_andrew) > 0)
    {
    gui_text_show(ERROR, " unknown scftyp ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  model->gamess.scf_type = bad_andrew;
  return(0);
  }

else if (g_ascii_strncasecmp(keyword, GMS_MAXIT_TXT, len=strlen(GMS_MAXIT_TXT)) == 0)
  model->gamess.maxit = (gint) (str_to_float(&keyword[len]));
else if (g_ascii_strncasecmp(keyword, GMS_TOTAL_Q_TXT, len=strlen(GMS_TOTAL_Q_TXT)) == 0)
  model->gamess.total_charge = (gint) (str_to_float(&keyword[len]));
else if (g_ascii_strncasecmp(keyword, GMS_MULT_TXT, len=strlen(GMS_MULT_TXT)) == 0)
  model->gamess.multiplicity = (gint) (str_to_float(&keyword[len]));
else if (g_ascii_strncasecmp(keyword, GMS_WIDE_OUTPUT_TXT, len=strlen(GMS_WIDE_OUTPUT_TXT)) == 0)
  model->gamess.wide_output = (((gint) (str_to_float(&keyword[len]))) == 6);
else if (g_ascii_strncasecmp(keyword, GMS_COORD_TXT, strlen(GMS_COORD_TXT)) == 0)//len UNUSED
  ; /* TODO handle different coordinate types */
else 
  {
  gui_text_show(ERROR, " unknown keyword ");
  gui_text_show(ERROR, keyword);
  gui_text_show(ERROR, " ");
  return(1);
  }
return(0);
}

gint get_system(gchar *keyword, struct model_pak *model)
{
int len;

if (g_ascii_strncasecmp(keyword, GMS_TIMLIM_TXT, len=strlen(GMS_TIMLIM_TXT)) == 0)
  model->gamess.time_limit = (gint) (str_to_float(&keyword[len]));
else if (g_ascii_strncasecmp(keyword, GMS_MWORDS_TXT, len=strlen(GMS_MWORDS_TXT)) == 0)
  model->gamess.mwords = (gint) (str_to_float(&keyword[len]));
else 
  {
  gui_text_show(ERROR, " unknown keyword ");
  gui_text_show(ERROR, keyword);
  gui_text_show(ERROR, " ");
  return(1);
  }
return(0);
}

gint get_statpt(gchar *keyword, struct model_pak *model)
{
int len, tmp;

if (g_ascii_strncasecmp(GMS_METHOD_TXT, keyword, len=strlen(GMS_METHOD_TXT)) == 0)
  {
  tmp = model->gamess.opt_type;
  if (read_keyword(&keyword[len], method_types, &tmp) > 0)
    {
    gui_text_show(ERROR, " unknown method ");
    gui_text_show(ERROR, &keyword[len]);
    gui_text_show(ERROR, " ");
    return(1);
    }
  return(0);
  }
else if (g_ascii_strncasecmp(keyword, GMS_NSTEP_TXT, len=strlen(GMS_NSTEP_TXT)) == 0)
  model->gamess.nstep = str_to_float(&keyword[len]);
else
  {
  gui_text_show(ERROR, " unknown keyword ");
  gui_text_show(ERROR, keyword);
  gui_text_show(ERROR, " ");
  return(1);
  }
return(0);
}


gint get_next_group(FILE *fp, struct model_pak *model, gint *have_basis)
{
gchar *line, *keyword;
gint ret;
GSList *keywords = NULL, *list;

line = file_read_line(fp);
if (!line)
  return(FALSE);

/* TODO not a valid keyword so for the moment skip but could store */    
if (g_ascii_strncasecmp(line, " $", 2) != 0)
  return(TRUE); 

if (g_ascii_strncasecmp(line, " $data", 6) == 0)
  {
  ret = get_data(fp, model, *have_basis);
  }
else if (g_ascii_strncasecmp(line, " $basis", 7) == 0)
  {
  *have_basis = TRUE;
  keywords = get_gamess_keywords(fp, line+7, &ret);
  for (list=keywords ; list ; list=g_slist_next(list))
    {
    keyword = (gchar *) list->data;
    ret = get_basis(keyword, model);  
    }
  }
else if (g_ascii_strncasecmp(line, " $contrl", 8) == 0)
  {
  keywords = get_gamess_keywords(fp, line+8, &ret);
  for (list=keywords; list ; list=g_slist_next(list))
    {
    keyword = (gchar *) list->data;
    ret = get_control(keyword, model);  
    }
  }
else if (g_ascii_strncasecmp(line, " $system", 8) == 0)
  {
  keywords = get_gamess_keywords(fp, line+8, &ret);
  for (list=keywords; list ; list=g_slist_next(list))
    {
    keyword = (gchar *) list->data;
    ret = get_system(keyword, model);  
    }
  }
else if (g_ascii_strncasecmp(line, " $statpt", 8) == 0)
  {
  keywords = get_gamess_keywords(fp, line+8, &ret);
  for (list=keywords; list ; list=g_slist_next(list))
    {
    keyword = (gchar *) list->data;
    ret = get_statpt(keyword, model);  
    }
  }
else if (g_ascii_strncasecmp(line, " $dft", 5) == 0)
  {
  model->gamess.dft = TRUE;
/* TODO - process functional */
  }
else
  {
  /* TODO - Unknown keyword, just pass through */
  }
free_slist(keywords);
g_free(line);
return(TRUE);
}

gint read_gms(gchar *filename, struct model_pak *model)
{
FILE *fp;
gchar line[LINELEN], *name;
gint have_basis = FALSE;

fp = fopen(filename, "rt");
if (!fp)
  {
  sprintf(line, "Unable to open file %s\n", filename);
  gui_text_show(ERROR, line);
  return(1);
  }
else
  {
  name = g_path_get_basename(filename);
  sprintf(line, "Opening %s: \n", name);
  g_free(name);
  gui_text_show(STANDARD, line);
  }

while (get_next_group(fp, model, &have_basis));

gui_text_show(STANDARD, "\n");

strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model_prep(model);
return(0);
}

/*******************************************/
/* read single GAMESS output configuration */
/*******************************************/
#define DEBUG_READ_GMS_OUT 1
gint read_gms_out_block(FILE *fp, struct model_pak *model, gint num_skip, gint bohr)
{
gint i, num_tokens;
gchar **buff, line[LINELEN];
GString *title, *energy_string, *grad_string;
GSList *clist;
struct core_pak *core;

clist = model->cores;

/* ignore first num_skip lines */
for (i=0 ; i<num_skip; i++)
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
    return(11);
    }

model->construct_pbc = FALSE;
model->fractional = FALSE;

/* get 1st line of coords */
if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
    return(11);
    }
buff = tokenize(line, &num_tokens);

while (num_tokens > 4)
  {
  if (clist)
    {
    core = clist->data;
    clist = g_slist_next(clist);
    }
  else
    {
    core = new_core(elements[(int) str_to_float(*(buff+1))].symbol, model);
    model->cores = g_slist_append(model->cores, core);
    }

  if (bohr)
    {
    core->x[0] = BOHR_TO_ANGS*str_to_float(*(buff+2));
    core->x[1] = BOHR_TO_ANGS*str_to_float(*(buff+3));
    core->x[2] = BOHR_TO_ANGS*str_to_float(*(buff+4));
    }
  else
    {
    core->x[0] = str_to_float(*(buff+2));
    core->x[1] = str_to_float(*(buff+3));
    core->x[2] = str_to_float(*(buff+4));
    }

/* get next line */
  g_strfreev(buff);
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading coordinates\n");
    return(11);
    }
  buff = tokenize(line, &num_tokens);
  }
g_strfreev(buff);
  
/* search for energy */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, " FINAL", 6) == 0)
    {
    buff = tokenize(line, &num_tokens);
    if (g_ascii_strncasecmp(*(buff+1), "ENERGY", 6) == 0)
      model->gamess.energy = str_to_float(*(buff+3));
    else
      model->gamess.energy = str_to_float(*(buff+4));
    model->gamess.have_energy = TRUE;
    g_strfreev(buff);
    break;
    }
  }

/* update for MP? */
if (model->gamess.MP_level > 0)
  {
  while (!fgetline(fp, line))
    {
    if (g_strrstr(line ,"E(MP2)") != NULL)
      {
      buff = tokenize(line, &num_tokens);
      model->gamess.energy = str_to_float(*(buff+1));
      model->gamess.have_energy = TRUE;
      g_strfreev(buff);
      break;
      }
    }
  }

/* search for gradient and read any properties */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, " NET CHARGES:", 13) == 0)
    {
    clist = model->cores;
    /* skip forward four lines */
    for (i=0 ; i<4; i++)
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading fitted charges\n");
        return(11);
        }
    while (clist != NULL)
      {
      buff = tokenize(line, &num_tokens);
      core = clist->data;
      core->lookup_charge = FALSE;
      core->charge = str_to_float(*(buff+1));
      g_strfreev(buff);
      clist = g_slist_next(clist);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading fitted charges\n");
        return(11);
        }
      }
    }
  if (g_ascii_strncasecmp(line, "          MAXIMUM GRADIENT", 26) == 0)
    {
    buff = tokenize(line, &num_tokens);
    model->gamess.max_grad = str_to_float(*(buff+3));
    model->gamess.have_max_grad = TRUE;
    model->gamess.rms_grad = str_to_float(*(buff+7));
    model->gamess.have_rms_grad = TRUE;
    g_strfreev(buff);
    /* check next line to see if converged */
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading equilibrium status\n");
      return(11);
      }
    if (g_ascii_strncasecmp(line, "1     ***** EQUILIBRIUM GEOMETRY LOCATED *****", 46) == 0)
      model->gamess.converged = TRUE;
    break;
    }
  }
g_free(model->title);
title = g_string_new("");
if (model->gamess.have_energy)
  {
  energy_string = g_string_new("");
  g_string_append_printf(energy_string, "%.5f H", model->gamess.energy);
  property_add_ranked(3, "Energy", energy_string->str, model);
  g_string_free(energy_string, TRUE); 
  g_string_append_printf(title, "E");
  if (model->gamess.MP_level > 0)
    g_string_append_printf(title, "(MP%d)", model->gamess.MP_level);
  g_string_append_printf(title, " = %.5f H", model->gamess.energy);
  }
if (model->gamess.have_rms_grad)
  {
  grad_string = g_string_new("");
  g_string_append_printf(grad_string, "%.4f H/B", model->gamess.rms_grad);
  property_add_ranked(4, "RMS Gradient", grad_string->str, model);
  g_string_free(grad_string, TRUE); 
  g_string_append_printf(title, ", grad = %.5f", model->gamess.rms_grad);
  }
if (model->gamess.have_max_grad)
  {
  grad_string = g_string_new("");
  g_string_append_printf(grad_string, "%.4f H/B", model->gamess.max_grad);
  property_add_ranked(4, "Maximum Gradient", grad_string->str, model);
  g_string_free(grad_string, TRUE); 
  }
model->title = g_strdup(title->str);
g_string_free(title, TRUE);

return(0);
}

/*******************************/
/* GAMESS output frame reading */
/*******************************/
gint read_gms_out_frame(FILE *fp, struct model_pak *model)
{
/* replace all data */
return(read_gms_out_block(fp, model, 2, FALSE));
}

/********************************/
/* Read in a GAMESS output file */
/********************************/
gint read_gms_out(gchar *filename, struct model_pak *model)
{
gint flag, frame, num_tokens, len, i, index, bad_andrew;
gchar **buff, line[LINELEN], *keyword, *option;
GString *property_string;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  {
  sprintf(line, "Unable to open file %s\n", filename);
  gui_text_show(ERROR, line);
  return(1);
  }

model->periodic = 0;
flag=frame=0;

/* read in BASIS OPTIONS */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "     BASIS OPTIONS", 18) == 0)
    {
    /* skip line */
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading basis options\n");
      return(2);
      }
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading basis options\n");
      return(2);
      }
    /* get first line of options i.e. basis set */
    buff = tokenize(line, &num_tokens); 
    /* GBASIS=STO          IGAUSS=       3      POLAR=NONE */
    keyword = *(buff+0);
    if (g_ascii_strncasecmp(keyword, GMS_BASIS_TXT, len = strlen(GMS_BASIS_TXT)) == 0)
      {
      if (read_keyword(&keyword[len], basis_types, &bad_andrew) > 0)
        {
        sprintf(line, "invalid basis %s\n", &keyword[len]);
        gui_text_show(ERROR, line);
        return(3);
        }
      model->gamess.basis = bad_andrew;
      }
    model->gamess.ngauss = (gint) str_to_float(*(buff+2));
    g_strfreev(buff);
    
    /* get 2nd line of options i.e. NDFUNC and DIFFSP */
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading basis options\n");
      return(2);
      }
    buff = tokenize(line, &num_tokens);
    /* NDFUNC=       0     DIFFSP=       F */
    model->gamess.num_d = str_to_float(*(buff+1));
    if (g_ascii_strncasecmp(*(buff+3), "F", 1) == 0)
      model->gamess.have_heavy_diffuse = FALSE;
    else
      model->gamess.have_heavy_diffuse = TRUE;
    g_strfreev(buff);
      
    /* get 3rd line of options i.e. MULT and ICHARG */
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading basis options\n");
      return(2);
      }
    buff = tokenize(line, &num_tokens);
    /* NPFUNC=       0      DIFFS=       F */
    model->gamess.num_p = (gint) str_to_float(*(buff+1));
    if (g_ascii_strncasecmp(*(buff+3), "F", 1) == 0)
      model->gamess.have_hydrogen_diffuse = FALSE;
    else
      model->gamess.have_hydrogen_diffuse = TRUE;
    g_strfreev(buff);
      
    /* TODO f functions */
	flag++;
	break;
	}
  }

if (!flag)
  {
   /* no basis present so set to user defined and rewind file */
   model->gamess.basis = GMS_USER;
   rewind(fp);
  }
flag=0;

/* read in RUN TITLE */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "     RUN TITLE", 14) == 0)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading title\n");
      return(2);
      }
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading title\n");
      return(2);
      }
    model->gamess.title = g_strdup(g_strstrip(line));
    flag++;
    break;
    }
  }

if (!flag)
  {
   gui_text_show(ERROR, "RUN TITLE not found\n");
   return(2);
  }
flag=0;

/* read in $CONTRL OPTIONS */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "     $CONTRL OPTIONS", 20) == 0)
    {
    flag++;
    if (fgetline(fp, line))
      /* skip line of dashes */
      {
      gui_text_show(ERROR, "unexpected end of file reading contrl options\n");
      return(3);
      }
    while (TRUE)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading contrl options\n");
        return(3);
        }
      /* is the line the blank line signalling end of control options? */
      if (strlen(g_strchug(line)) == 0)
        break;
      /* break up line into option pairs */
      /* each pair takes 15 characters with 5 characters between them */
      /* note that we have already removed the single space in front of the lines with the g_strchug */
      index = 0;
      while (index+15 <= strlen(line))
        {
        option = g_strndup(line+index, 15);
        /* split into pair */
        buff = g_strsplit(option, "=", 2);
        g_free(option);
        /* remove whitespace */
        g_strstrip(buff[0]);
        g_strstrip(buff[1]);
        /* the compare strings end in = which we have stripped off so compare on strlen-1 */
        if (g_ascii_strncasecmp(buff[0], GMS_SCFTYPE_TXT, strlen(GMS_SCFTYPE_TXT) - 1) == 0)
          {
          if (read_keyword(buff[1], scf_types, &bad_andrew) > 0)
            {
            sprintf(line, "invalid scf type %s\n", buff[1]);
            gui_text_show(ERROR, line);
            return(3);
            }
          model->gamess.scf_type = bad_andrew;
          }
        else if (g_ascii_strncasecmp(buff[0], GMS_RUNTYPE_TXT, strlen(GMS_RUNTYPE_TXT) - 1) == 0)
          {
          if (read_keyword(buff[1], run_types, &bad_andrew) > 0)
            {
            sprintf(line, "invalid run type %s\n", buff[1]);
            gui_text_show(ERROR, line);
            return(3);
            }
          model->gamess.run_type = bad_andrew;
            property_string = g_string_new("");
            i=0;
            while (run_types[i].label)
              {
              if (model->gamess.run_type == run_types[i].id)
                g_string_append_printf(property_string, "%s", run_types[i].label);
              i++;
              }
            property_string->str[0] = g_ascii_toupper(property_string->str[0]);
            property_add_ranked(2, "Calculation", property_string->str, model);
            g_string_free(property_string, TRUE); 
          }
        else if (g_ascii_strncasecmp(buff[0], GMS_EXETYPE_TXT, strlen(GMS_EXETYPE_TXT) - 1) == 0)
          {
          if (read_keyword(buff[1], exe_types, &bad_andrew) > 0)
            {
            sprintf(line, "invalid execution type %s\n", buff[1]);
            gui_text_show(ERROR, line);
            return(3);
            }
          model->gamess.exe_type = bad_andrew;
          }
        else if (g_ascii_strncasecmp(buff[0], GMS_MPLEVEL_TXT, strlen(GMS_MPLEVEL_TXT) - 1) == 0)
        	model->gamess.MP_level = (gint) str_to_float(buff[1]);
        else if (g_ascii_strncasecmp(buff[0], GMS_CITYP_TXT, strlen(GMS_CITYP_TXT) - 1) == 0)
          if (g_ascii_strncasecmp(buff[1], "none", 4) == 0)
            model->gamess.have_CI = FALSE;
          else
            model->gamess.have_CI = TRUE;
        else if (g_ascii_strncasecmp(buff[0], GMS_CCTYP_TXT, strlen(GMS_CCTYP_TXT) - 1) == 0)
          if (g_ascii_strncasecmp(buff[1], "none", 4) == 0)
            model->gamess.have_CC = FALSE;
          else
            model->gamess.have_CC = TRUE;
        else if (g_ascii_strncasecmp(buff[0], GMS_TOTAL_Q_TXT, strlen(GMS_TOTAL_Q_TXT) - 1) == 0)
          model->gamess.total_charge = (gint) str_to_float(buff[1]);
        else if (g_ascii_strncasecmp(buff[0], GMS_MULT_TXT, strlen(GMS_MULT_TXT) - 1) == 0)
          model->gamess.multiplicity = (gint) str_to_float(buff[1]);
        else if (g_ascii_strncasecmp(buff[0], GMS_MAXIT_TXT, strlen(GMS_MAXIT_TXT) - 1) == 0)
          model->gamess.maxit = ((gint) str_to_float(buff[1]));
        else if (g_ascii_strncasecmp(buff[0], GMS_WIDE_OUTPUT_TXT, strlen(GMS_WIDE_OUTPUT_TXT) - 1) == 0)   
         model->gamess.wide_output = ((gint) str_to_float(buff[1]) == 6);

        g_strfreev(buff);
        index += 20;
        }
      }
      break;
    }
  }

if (!flag)
  {
/* don't return... model_prep() needs to be called to avoid crashing */
   gui_text_show(WARNING, "$CONTRL OPTIONS not found\n");
  }
flag=0;

/* read in $SYSTEM OPTIONS */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "     $SYSTEM OPTIONS", 20) == 0)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading system options\n");
      return(4);
      }
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading system options\n");
      return(4);
      }
    buff = tokenize(line, &num_tokens);
    model->gamess.mwords = (gint) (str_to_float(*(buff+2))/1000000);
    g_strfreev(buff);
    
    for (i=0; i<4; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading system options\n");
        return(4);
        }
      }
    buff = tokenize(line, &num_tokens);
    model->gamess.time_limit = (gint) (str_to_float(*(buff+1))/60.0);
    g_strfreev(buff);
    flag++;
    break;
    }
  }

if (!flag)
  {
/* don't return... model_prep() needs to be called to avoid crashing */
   gui_text_show(WARNING, "$SYSTEM OPTIONS not found\n");
  }
flag=0;

/* anything else to find ? */
while (!fgetline(fp, line))
  {
  if (g_ascii_strncasecmp(line, "                         GRADIENT OF THE ENERGY", 47) == 0)
    {
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading gradient\n");
      return(5);
      }
    while (g_ascii_strncasecmp(line, "                   MAXIMUM GRADIENT", 35) != 0) 
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading gradient\n");
        return(5);
        }
      }
      buff = tokenize(line, &num_tokens);
      model->gamess.max_grad = str_to_float(*(buff+3));
      model->gamess.have_max_grad = TRUE;
      g_strfreev(buff);
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading gradient\n");
        return(5);
        }
      buff = tokenize(line, &num_tokens);
      model->gamess.rms_grad = str_to_float(*(buff+3));
      model->gamess.have_rms_grad = TRUE;
      g_strfreev(buff);
    }
  }

rewind(fp);

/* Read the input coordinates - single frame has different format to multiframe */
if (model->gamess.run_type < GMS_OPTIMIZE) { /* is it a single frame job? */
  while (!fgetline(fp, line))
    {
    if (g_ascii_strncasecmp(line, " ATOM      ATOMIC                      COORDINATES (BOHR)", 57) == 0)
      {
      read_gms_out_block(fp, model, 1, TRUE);
      flag++;
      break;
      }
    }
  }

else
  {
  /* get optimisation parameters */
  while (!fgetline(fp, line))
    {
    if (g_ascii_strncasecmp(line, "          STATIONARY POINT LOCATION RUN", 39) == 0)
      {
      for (i=0; i<7; i++)
        {
        if (fgetline(fp, line))
          {
          gui_text_show(ERROR, "unexpected end of file reading optimizer options\n");
          return(5);
          }
        }
      if (model->gamess.exe_type == GMS_CHECK)
        if (fgetline(fp, line))
          {
          gui_text_show(ERROR, "unexpected end of file reading optimizer options\n");
          return(5);
          }
      buff = tokenize(line, &num_tokens);
      if (read_keyword(&(*(buff+1))[1], method_types, &bad_andrew) > 0)
        {
        sprintf(line, "invalid method %s\n",&(*(buff+1))[1]);
        gui_text_show(ERROR, line);
        return(5);
        }
      model->gamess.opt_type = bad_andrew;
      g_strfreev(buff);
      flag++;
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading optimizer options\n");
        return(5);
        }
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading optimizer options\n");
        return(5);
        }
      buff = tokenize(line, &num_tokens);
      model->gamess.nstep = str_to_float(*(buff+2));
      g_strfreev(buff);
      flag++;
      break;
      }
    }
  if (!flag)
    {
    gui_text_show(ERROR, "optimizer options not found\n");
    return(5);
    }
  /* Are there any coordinates from a minimisation? */
  flag=0;
  while (!fgetline(fp, line) && !model->gamess.converged)
    {
    /* coordinates */
    if (g_ascii_strncasecmp(line, " COORDINATES OF ALL ATOMS ARE", 29) == 0)
      {
      /* go through all frames to count them */
      add_frame_offset(fp, model);
      read_gms_out_block(fp, model, 2, FALSE);
      flag++;
      frame++;
      }
    }
  }

/*  property_string = g_string_new("");
  g_string_append_printf(property_string, "%.0f", model->gamess.total_charge);
  property_add_ranked(2, "Total Charge", property_string->str, model);
  g_string_free(property_string, TRUE); 
  property_string = g_string_new("");
  g_string_append_printf(property_string, "%.0f", model->gamess.multiplicity);
  property_add_ranked(3, "Multiplicity", property_string->str, model);
  g_string_free(property_string, TRUE); */

  property_string = g_string_new("");
  i=0;
  while (basis_sets[i].label)
    {
    if ((basis_sets[i].basis == model->gamess.basis) && (basis_sets[i].ngauss == model->gamess.ngauss))
      g_string_append_printf(property_string, "%s", basis_sets[i].label);
    i++;
    }
  property_add_ranked(7, "Basis", property_string->str, model);
  g_string_free(property_string, TRUE); 

/* done */

if (flag)
  {
  /* set frame if don't want last? */
  strcpy(model->filename, filename);
  g_free(model->basename);
  model->basename = parse_strip(filename);

  model->num_frames = model->cur_frame = frame;
  model->cur_frame--;
  
  model_prep(model);
  }
else
  return(2);

return(0);
}

