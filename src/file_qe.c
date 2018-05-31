/*
Copyright (C) 2015 by Andrew Lloyd Rohl

a.rohl@curtin.edu.au

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

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "parse.h"
#include "matrix.h"
#include "model.h"
#include "interface.h"
#include "math.h"

#define DEBUG_READ_QE_OUT 0
#define DEBUG_READ_QE_IN 0

typedef enum {
  QE_CELL_ALAT, QE_CELL_BOHR /*, QE_CELL_ANGSTROM*/
} QECellUnits;

typedef enum {
  QE_COORD_ALAT, QE_COORD_FRAC, QE_COORD_ANGSTROM, QE_COORD_BOHR
} QECoordUnits;

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* global variables */

gdouble alat;


/************************************/
/* Quantum ESPRESSO compare species */
/************************************/
gint comp_species(gpointer pa, gpointer pb)
{
  const struct qe_species_pak *a = pa, *b = pb;
  
  /* Compare by label */
  return strcmp(a->label, b->label);
}

/********************************/
/* Quantum ESPRESSO file saving */
/********************************/
gint write_qe(gchar *filename, struct model_pak *model)
{
gint i;
gdouble x[3];
GSList *list;
struct core_pak *core;
struct qe_species_pak *species;
gchar *element;
FILE *fp;
  
/* checks */
g_return_val_if_fail(model != NULL, 1);
g_return_val_if_fail(filename != NULL, 2);
  
/* open the file */
fp = fopen(filename,"wt");
if (!fp)
  return(3);

if (model->periodic != 3)
  {
  gui_text_show(ERROR, "Your model must have 3D periodicity\n");
  return(2);
  }

if (model->qe.species == NULL)
  {
  /* No species data - assume that species are simply the elements */

  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = (struct core_pak *) list->data;
    if (core->status & DELETED)
      continue;
    element = elements[core->atom_code].symbol;
    species = g_malloc(sizeof(struct qe_species_pak));
    species->label = g_strdup(element);
    /* if this species doesn't exist in list, add */
    if (g_slist_find_custom(model->qe.species, species, (GCompareFunc)comp_species) == NULL)
      {
      species->mass = elements[core->atom_code].weight;
      species->pseudo = g_strdup_printf("%s.UPF", element);
#ifdef DEBUG_READ_QE_OUT
      printf("adding %s\n", species->label);
#endif
      model->qe.species = g_slist_append(model->qe.species, species);
      }
    }
  }
  
g_free(model->basename);
model->basename = parse_strip(filename);
  
fprintf(fp, "&control\n");
fprintf(fp, "  title = '%s',\n", model->title);
fprintf(fp, "  prefix = '%s',\n", model->basename);
fprintf(fp, "  pseudo_dir='.',\n");
fprintf(fp, "/\n");
  
fprintf(fp, "\n&system\n");
fprintf(fp, "  ibrav = 0,\n");
fprintf(fp, "  nat = %d,\n", model->num_atoms);
fprintf(fp, "  ntyp = %d,\n", g_slist_length(model->qe.species));
fprintf(fp, "  ecutwfc = %11.4f,\n", model->qe.ecutwfc);
fprintf(fp, "/\n");
  
fprintf(fp, "\n&electrons\n");
fprintf(fp, "/\n");

fprintf(fp, "\n&ions\n");
fprintf(fp, "/\n");

fprintf(fp, "\n&cell\n");
fprintf(fp, "/\n");

fprintf(fp, "\nATOMIC_SPECIES\n");
for (list=model->qe.species; list; list=g_slist_next(list))
  {
  species = (struct qe_species_pak *) list->data;
  fprintf(fp, "%-4s %10.6f %s\n", species->label, species->mass, species->pseudo);
  }
  
fprintf(fp, "\nCELL_PARAMETERS (angstrom)\n");
for (i=0; i<model->periodic; i++)
  fprintf(fp, "%14.9f  %14.9f  %14.9f\n", model->latmat[0+i], model->latmat[3+i], model->latmat[6+i]);

fprintf(fp, "\nATOMIC_POSITIONS (angstrom)\n");
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = (struct core_pak *) list->data;
  if (core->status & DELETED)
    continue;
    
  /* everything is cartesian after latmat mult */
  ARR3SET(x, core->x);
  vecmat(model->latmat, x);
  fprintf(fp, "%-7s    %14.9f  %14.9f  %14.9f\n", elements[core->atom_code].symbol, x[0], x[1], x[2]);
  }
fprintf(fp, "\n");
fclose(fp);
return(0);
}

//gchar **strip_empty_strings(gchar **str_array)
//{
//gchar **new_str_array = NULL;
//gint i=0, j=0;
//
//while (str_array[i] != NULL)
//  {
//  if (strlen(str_array[i]) != 0)
//    {
//    new_str_array = g_realloc(new_str_array, (j + 2)*sizeof(char *)); /* add one for NULL pointer */
//    new_str_array[j] = str_array[i];
//    j++;
//    }
//  else
//    g_free(str_array[i]);
//  i++;
//  }
//  new_str_array[j] = NULL;
//  g_free(str_array);
//  return(new_str_array);
//}
//
//
//gchar **read_qe_keyword(gchar *string)
//{
//  gint i;
///* 
//Options to keywords may be specified in any of the following forms:
//
//            keyword = option  
//            keyword(option) 
//            keyword=(option1, option2, ...)  
//            keyword(option1, option2, ...) 
//*/
//
///* Spaces, tabs, commas, or forward slashes can be used in any combination to separate items within a line. Multiple spaces are treated as a single delimiter. */
//
//for (i=0; i<strlen(string); i++)
//  switch (string[i])
//    {
//    case '=':
//    case '(':
//    case ')':
//    case '\t':
//    case ',': string[i] = ' ';
//    }
//return(strip_empty_strings(g_strsplit(string, " ", 0)));
//}

/**********************************/
/* process an ATOMIC_SPECIES line */
/**********************************/
//void qe_parse_species(const gchar *line, GSList **list)
//{
//  gint num_tokens;
//  gchar **buff;
//  struct qe_species_pak *species;
//  
//  g_assert(line != NULL);
//  
//  buff = tokenize(line, &num_tokens);
//  
//  if (num_tokens > 2)
//    {
//    species = g_malloc(sizeof(struct qe_species_pak));
//    
//    species->label = g_strdup(*(buff+2));
//    *list = g_slist_append(*list, species);
//    }
//  
//  g_strfreev(buff);
//}


/*****************************************************/
/* read single Quantum ESPRESSO output configuration */
/*****************************************************/
gint read_qe_out_block(FILE *fp, struct model_pak *model)
{
gint i;
gint num_tokens;
gchar **buff, line[LINELEN], *sub_line, *text;
GSList *clist, *found;
GString *grad_string, *stress_string;
struct core_pak *core;
struct qe_species_pak *species, *found_species;
  
/* TODO: check for final frame... */
  
clist = model->cores;

model->periodic = 3;
model->construct_pbc = TRUE;

/* read to end of iteration */
while (TRUE)
  {
  if (fgetline(fp, line))
    {
    gui_text_show(ERROR, "unexpected end of file reading to end of iteration\n");
    return(2);
    }
  
  /* read in cell vectors of first/last frame) */
  if (g_strrstr(line, "crystal axes:") != NULL)
    {
    /* NB: gdis wants transposed matrix */
   for (i=0; i<3; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading cell vectors\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      model->latmat[0+i] = str_to_float(*(buff+3))*alat*BOHR_TO_ANGS;
      model->latmat[3+i] = str_to_float(*(buff+4))*alat*BOHR_TO_ANGS;
      model->latmat[6+i] = str_to_float(*(buff+5))*alat*BOHR_TO_ANGS;
      g_strfreev(buff);
      }
    }
  
  /* read in pseudos - note that in a conp optimisation, this data appears in first and last frame */
  if ((g_strrstr(line, "PseudoPot.") != NULL) && !model->qe.have_species)
    {
    species = g_malloc(sizeof(struct qe_species_pak));
    sub_line = g_strrstr(line, "for ");
    buff = tokenize((gchar *)(sub_line), &num_tokens);
    species->label = g_strdup(buff[1]);
    g_strfreev(buff);
    if (fgetline(fp, line))
      {
      gui_text_show(ERROR, "unexpected end of file reading name of pseudopotential file\n");
      return(2);
      }
    buff = tokenize((gchar *)(line), &num_tokens);
    species->pseudo = g_strdup(g_path_get_basename(buff[0]));
    model->qe.species = g_slist_append(model->qe.species, species);
    g_strfreev(buff);
    }
  
  /* read in atomic data */
  if (g_strrstr(line, "atomic species") != NULL)
    {
    while (!fgetline(fp, line))
      {
      species = g_malloc(sizeof(struct qe_species_pak));
      buff = tokenize((gchar *)(line), &num_tokens);
      if (num_tokens == 0)
        {
        model->qe.have_species = TRUE;
        break;
        }
      else if (num_tokens == 6)
        {
        species = g_malloc(sizeof(struct qe_species_pak));
        species->label = g_strdup(buff[0]);
        if ((found = g_slist_find_custom(model->qe.species, species, (GCompareFunc)comp_species)) == NULL)
          {
          gui_text_show(ERROR, "unexpected end of file reading atomic species data\n");
          g_free(species);
          return(2);
          }
        else
          {
          found_species = (struct qe_species_pak *) found->data;
          found_species->mass = str_to_float(*(buff+2));
          g_free(species);
          }
        }
      else
        {
        gui_text_show(ERROR, "unexpected number of tokens reading atomic species data (6 expected)\n");
        return(2);
        }
      g_strfreev(buff);
      }
    }

  /* read in cell vectors of other frames */
  if (g_strrstr(line, "CELL_PARAMETERS") != NULL)
    {
    QECellUnits cell_units;
    gdouble multiplier;
    
    buff = tokenize(line, &num_tokens);
    if (num_tokens == 1)
      cell_units = QE_CELL_BOHR;
    else
      {
      strip_extra(*(buff+1));
      if (g_ascii_strncasecmp(*(buff+1), "alat", 4) == 0)
        cell_units = QE_CELL_ALAT;
      else if (g_ascii_strncasecmp(*(buff+1), "bohr", 4) == 0)
        cell_units = QE_CELL_BOHR;
  //    else if (g_ascii_strncasecmp(*(buff+1), "angstrom", 8) == 0)
  //      cell_units = QE_CELL_ANGSTROM;
      else
        {
        gchar *msg;
        
        msg = g_strdup_printf("%s: Unknown units for cell parameters\n", *(buff+1));
        gui_text_show(ERROR, msg);
        g_free(msg);
        return(2);
        }
      }
    g_strfreev(buff);
    switch (cell_units)
      {
        case QE_CELL_BOHR:
          multiplier = BOHR_TO_ANGS;
        break;
        case QE_CELL_ALAT:
          multiplier = alat*BOHR_TO_ANGS;
        break;
      }
    /* NB: gdis wants transposed matrix */
    for (i=0; i<3; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file reading cell vectors\n");
        return(2);
        }
      buff = tokenize(line, &num_tokens);
      model->latmat[0+i] = str_to_float(*(buff+0))*multiplier;
      model->latmat[3+i] = str_to_float(*(buff+1))*multiplier;
      model->latmat[6+i] = str_to_float(*(buff+2))*multiplier;
      g_strfreev(buff);
      }
    }

  /* read coordinates of first/last frame */
  if (g_strrstr(line, "Cartesian axes") != NULL)
    {
    model->fractional = FALSE;
    /* skip 2 lines */
    for (i=0; i<2; i++)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file skipping lines in reading coordinates\n");
        return(2);
        }
      }
    while (!fgetline(fp, line))
      {
      buff = tokenize(line, &num_tokens);
      if (num_tokens == 10)
        {
        if (clist)
          {
          core = (struct core_pak *) clist->data;
          clist = g_slist_next(clist);
          }
        else
          {
          core = new_core(*(buff+1), model);
          model->cores = g_slist_append(model->cores, core);
          }
        core->x[0] = str_to_float(*(buff+6))*alat*BOHR_TO_ANGS;
        core->x[1] = str_to_float(*(buff+7))*alat*BOHR_TO_ANGS;
        core->x[2] = str_to_float(*(buff+8))*alat*BOHR_TO_ANGS;
#if DEBUG_READ_QE_OUT
        printf("new coords %f %f %f\n", core->x[0], core->x[1], core->x[2]);
#endif
        g_strfreev(buff);
        }
      else
        {
        g_strfreev(buff);
        break;
        }
      }
    }
  
  /* read coordinates of other frames */
  if (g_strrstr(line, "ATOMIC_POSITIONS") != NULL)
    {
    QECoordUnits coord_units;
    gdouble multiplier;
    
    buff = tokenize(line, &num_tokens);
    if (num_tokens == 1)
      coord_units = QE_COORD_BOHR;
    else
      {
      strip_extra(*(buff+1));
      if (g_ascii_strncasecmp(*(buff+1), "alat", 4) == 0)
        coord_units = QE_COORD_ALAT;
      else if (g_ascii_strncasecmp(*(buff+1), "crystal", 7) == 0)
        coord_units = QE_COORD_FRAC;
      else if (g_ascii_strncasecmp(*(buff+1), "angstrom", 8) == 0)
        coord_units = QE_COORD_ANGSTROM;
      else
        {
        gchar *msg;
        
        msg = g_strdup_printf("%s: Unknown units for coordinates\n", *(buff+1));
        gui_text_show(ERROR, msg);
        g_free(msg);
        return(2);
        }
      }
    g_strfreev(buff);
    switch (coord_units)
      {
        case QE_COORD_ALAT:
          model->fractional = FALSE;
          multiplier = alat*BOHR_TO_ANGS;
        break;
        case QE_COORD_ANGSTROM:
          model->fractional = FALSE;
          multiplier = 1.0;
        break;
        case QE_COORD_BOHR:
          model->fractional = FALSE;
          multiplier = BOHR_TO_ANGS;
        break;
        case QE_COORD_FRAC:
          model->fractional = TRUE;
          multiplier = 1.0;
        break;
      }

    while (!fgetline(fp, line))
      {
      buff = tokenize(line, &num_tokens);
      if (num_tokens == 4)
        {
        if (clist)
          {
          core = (struct core_pak *) clist->data;
          clist = g_slist_next(clist);
          }
        else
          {
          core = new_core(*(buff+0), model);
          model->cores = g_slist_append(model->cores, core);
          }
        core->x[0] = str_to_float(*(buff+1))*multiplier;
        core->x[1] = str_to_float(*(buff+2))*multiplier;
        core->x[2] = str_to_float(*(buff+3))*multiplier;
#if DEBUG_READ_QE_OUT
        printf("new coords %f %f %f\n", core->x[0], core->x[1], core->x[2]);
#endif
        g_strfreev(buff);
        }
      else
        {
        g_strfreev(buff);
        break;
        }
      }
    }
  
  /* CPMD specific quantities */
  
  /* energy */
  if (g_strrstr(line, "!") != NULL)
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    g_strfreev(buff);
    }
  
  /* force */
  if (g_strrstr(line, "Total force") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    grad_string = g_string_new("");
    g_string_append_printf(grad_string, "%s %s", *(buff+3), "Ry/au");
    text = format_value_and_units(grad_string->str, 5);
    property_add_ranked(4, "Force", text, model);
    g_free(text);
    g_string_free(grad_string, TRUE);
    g_strfreev(buff);
    }
  
  /* stress */
  if (g_strrstr(line, "total   stress") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    stress_string = g_string_new("");
    g_string_append_printf(stress_string, "%s %s", *(buff+5), "kbar");
    text = format_value_and_units(stress_string->str, 2);
    property_add_ranked(5, "Stress", text, model);
    g_free(text);
    g_string_free(stress_string, TRUE);
    g_strfreev(buff);
    }
  
  /* Quantum ESPRESSO writes out the final coordinates twice - skip second set */
  if (g_strrstr(line, "Begin final coordinates") != NULL)
    {
    while (TRUE)
      {
      if (fgetline(fp, line))
        {
        gui_text_show(ERROR, "unexpected end of file skipping to end of final coordinates repeat\n");
        return(2);
        }
      if (g_strrstr(line, "End final coordinates") != NULL)
        break;
      }
    }

  /* CPMD specific quantities */
  
  /* energy */
  if (g_strrstr(line, "total energy") != NULL)
    {/* _BUG_ this will catch "The total energy is the sum of (...)" line as well*/
 if(g_strrstr(line, "total energy is the sum of") == NULL)/* FIX: exclude this line */
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 5);
    property_add_ranked(3, "Energy", text, model);
    g_free(text);
    g_strfreev(buff);
    }
 }

  /* end of first frame */
  if (g_strrstr(line, "init_run") != NULL)
    break;
  
  /* end of other frames */
  if (g_strrstr(line, "number of scf cycles    ") != NULL)
    break;
  
  /* end of CP MD frame */
  if (g_strrstr(line, "MLWF step") != NULL)
    break;

  /* end of optimisation */
  if (g_strrstr(line, "A final scf calculation at the relaxed structure") != NULL)
    break;
  }

#if DEBUG_READ_QE_OUT
  printf("--- end of frame ---\n");
#endif
return(0);
}

/*****************************************/
/* Quantum ESPRESSO output frame reading */
/*****************************************/
gint read_qe_out_frame(FILE *fp, struct model_pak *model)
{
    return(read_qe_out_block(fp, model));
}

/******************************************/
/* Read in a Quantum ESPRESSO output file */
/******************************************/
gint read_qe_out(gchar *filename, struct model_pak *model)
{
gint frame, num_tokens;
gchar **buff, *text, line[LINELEN];
FILE *fp;
gint last_frame;
  
fp = fopen(filename, "rt");
if (!fp)
  return(1);

last_frame = FALSE;
frame=0;


while (!fgetline(fp, line))
  {
  if (g_strrstr(line, "Program PWSCF") != NULL)
    {
    property_add_ranked(2, "Calculation", "PWSCF", model);
    }
  if (g_strrstr(line, "Program CP") != NULL)
    {
      property_add_ranked(2, "Calculation", "CPMD", model);
    }
  if (g_strrstr(line, "Exchange-correlation") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    property_add_ranked(7, "Ex-Corr", buff[2], model);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "EXX-fraction") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    property_add_ranked(8, "Ex-fraction", buff[2], model);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "kinetic-energy cutoff") != NULL)
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 2);
    property_add_ranked(6, "KE cutoff", text, model);
    g_free(text);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "Ecutwfc") != NULL)
    {
    buff = g_strsplit(line, "=", 2);
    text = format_value_and_units(*(buff+1), 2);
    property_add_ranked(6, "KE cutoff", text, model);
    g_free(text);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "number of k points") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    property_add_ranked(10, "K-points", buff[4], model);
    g_strfreev(buff);
    }
  if (g_strrstr(line, "lattice parameter (alat)") != NULL)
    {
    buff = tokenize((gchar *)(line), &num_tokens);
    alat = str_to_float(*(buff+4));
    g_strfreev(buff);
    }
  if ((g_strrstr(line, "number of bfgs") != NULL) || (g_strrstr(line, "celldm") != NULL) || (g_strrstr(line, "Physical Quantities at step") != NULL))
    {
		/* go through all frames to count them */
		add_frame_offset(fp, model);
		read_qe_out_block(fp, model);
		frame++;
    if (last_frame)
      break;
    }
  }

/* done */
strcpy(model->filename, filename);
g_free(model->basename);
model->basename = parse_strip(filename);
model->num_frames = model->cur_frame = frame;
model->cur_frame--;

model_prep(model);

return(0);
}


/*********************************************/
/* get a non trivial line from QE input file */
/*********************************************/
gint fgetqeline(FILE *fp, gchar *line)
{
  gint i, linlen;
  gboolean found_comment;
  
    /* get a line */
  if (fgets(line, LINELEN/2, fp) == NULL)
    return(1);
  linlen = strlen(line);
  /* only treated as data if not a comment and non empty */
  if (line[0] != '!' && linlen)
    {
    /* Anything following a ! is a comment so remove */
    found_comment = FALSE;
    for (i=0; i<strlen(line); i++)
      {
      if (line[i] == '!')
        found_comment = TRUE;
      if (found_comment)
        line[i] = ' ';
      /* strip quotes and commas */
      if (line[i] == ',' || line[i] == '\"' || line[i] == '\'')
        line[i] = ' ';
      }
    g_strstrip(line);
    }
  /* all clear */
  return(0);
}


/******************************************/
/* Read in a Quantum ESPRESSO input file */
/******************************************/
gint read_qe(gchar *filename, struct model_pak *model)
{
  gint i, num_tokens;
  gint nat=0, ibrav=0;
#ifdef UNUSED_BUT_SET
  gint ntyp=0;
#endif
  gdouble celldm[6] ={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  gchar **buff, *text, line[LINELEN];
  FILE *fp;
  GSList *clist;
  struct core_pak *core;

  fp = fopen(filename, "rt");
  if (!fp)
    return(1);
  
  clist = model->cores;
  
  model->periodic = 3;
  
  while (!fgetqeline(fp, line))
    {
    /* convert line to lower case */
    for (i=0; i<strlen(line); i++)
      {
      line[i] = g_ascii_tolower(line[i]);
      }

    if (g_strrstr(line, "&system") != NULL)
      {
      while (TRUE)
        {
        /* TODO: should case be changed? */
        /* or write own strrstr which copies strings, converts both to lowercase and then does the strstr */
        fgetqeline(fp, line);
        if (line[0] == '/')
          break;
        if (g_strrstr(line, "nat") != NULL)
          {
          buff = g_strsplit(line, "=", 2);
          nat = (gint) str_to_float(*(buff+1));
          g_strfreev(buff);
          }
        else if (g_strrstr(line, "ntyp") != NULL)
          {
          buff = g_strsplit(line, "=", 2);
#ifdef UNUSED_BUT_SET
          ntyp = (gint) str_to_float(*(buff+1));
#endif
          g_strfreev(buff);
          }
        else if (g_strrstr(line, "ecutwfc") != NULL)
          {
          buff = g_strsplit(line, "=", 2);
          model->qe.ecutwfc = str_to_float(*(buff+1));
          text = format_value_and_units(*(buff+1), 2);
          property_add_ranked(6, "KE cutoff", text, model);
          g_free(text);
          g_strfreev(buff);
          }
        else if (g_strrstr(line, "input_dft") != NULL)
          {
          buff = g_strsplit(line, "=", 2);
          model->qe.input_dft = g_strdup(*(buff+1));
          property_add_ranked(7, "Ex-Corr", model->qe.input_dft, model);
          g_strfreev(buff);
          }
        else if (g_strrstr(line, "celldm") != NULL)
          {
          gint index;
          
          index = line[7] - '1';
          buff = g_strsplit(line, "=", 2);
          celldm[index] = str_to_float(*(buff+1));
          g_strfreev(buff);
          }
        else if (g_strrstr(line, "ibrav") != NULL)
          {
          buff = g_strsplit(line, "=", 2);
          ibrav = (gint) str_to_float(*(buff+1));
          g_strfreev(buff);
          }
        }
      }
    if (g_strrstr(line, "atomic_positions") != NULL)
      {
      QECoordUnits coord_units;
      gdouble multiplier;
      
      if (nat < 1)
        {
        gui_text_show(ERROR, "number of atoms not defined or less than 1");
        return(2);
        }
        
      buff = tokenize(line, &num_tokens);
      if (num_tokens == 1)
        coord_units = QE_COORD_BOHR;
      else
        {
        strip_extra(*(buff+1));
        if (g_ascii_strncasecmp(*(buff+1), "alat", 4) == 0)
          coord_units = QE_COORD_ALAT;
        else if (g_ascii_strncasecmp(*(buff+1), "crystal", 7) == 0)
          coord_units = QE_COORD_FRAC;
        else if (g_ascii_strncasecmp(*(buff+1), "angstrom", 8) == 0)
          coord_units = QE_COORD_ANGSTROM;
        else if (g_ascii_strncasecmp(*(buff+1), "bohr", 8) == 0)
          coord_units = QE_COORD_BOHR;
        else
          {
          gchar *msg;
          
          msg = g_strdup_printf("%s: Unknown units for coordinates\n", *(buff+1));
          gui_text_show(ERROR, msg);
          g_free(msg);
          return(2);
          }
        }
      g_strfreev(buff);
      switch (coord_units)
        {
          case QE_COORD_ALAT:
          model->fractional = FALSE;
          multiplier = alat*BOHR_TO_ANGS;
          break;
          case QE_COORD_ANGSTROM:
          model->fractional = FALSE;
          multiplier = 1.0;
          break;
          case QE_COORD_BOHR:
          model->fractional = FALSE;
          multiplier = BOHR_TO_ANGS;
          break;
          case QE_COORD_FRAC:
          model->fractional = TRUE;
          multiplier = 1.0;
          break;
        }
      for (i=0; i< nat; i++)
        {
        fgetqeline(fp, line);
        buff = tokenize(line, &num_tokens);
        if (num_tokens == 4)
          {
          if (clist)
            {
            core = (struct core_pak *) clist->data;
            clist = g_slist_next(clist);
            }
          else
            {
            core = new_core(*(buff+0), model);
            model->cores = g_slist_append(model->cores, core);
            }
          core->x[0] = str_to_float(*(buff+1))*multiplier;
          core->x[1] = str_to_float(*(buff+2))*multiplier;
          core->x[2] = str_to_float(*(buff+3))*multiplier;
#if DEBUG_READ_QE_OUT
          printf("new coords %f %f %f\n", core->x[0], core->x[1], core->x[2]);
#endif
          g_strfreev(buff);
          }
        else
          {
          g_strfreev(buff);
          break;
          }
        }
      }
    if (g_strrstr(line, "cell_parameters") != NULL)
      {
      QECoordUnits cell_units;
      gdouble multiplier;

      buff = tokenize(line, &num_tokens);
      if (num_tokens == 1)
        cell_units = QE_COORD_BOHR;
      else
        {
        strip_extra(*(buff+1));
        if (g_ascii_strncasecmp(*(buff+1), "alat", 4) == 0)
          cell_units = QE_COORD_ALAT;
        else if (g_ascii_strncasecmp(*(buff+1), "angstrom", 8) == 0)
          cell_units = QE_COORD_ANGSTROM;
        else if (g_ascii_strncasecmp(*(buff+1), "bohr", 8) == 0)
          cell_units = QE_COORD_BOHR;
        else
          {
          gchar *msg;
          
          printf("%s: Unknown units for cell\n", *(buff+1));
          msg = g_strdup_printf("%s: Unknown units for cell\n", *(buff+1));
          gui_text_show(ERROR, msg);
          g_free(msg);
          return(2);
          }
        }
      g_strfreev(buff);
      switch (cell_units)
        {
          case QE_COORD_ALAT:
          multiplier = alat*BOHR_TO_ANGS;
          break;
          case QE_COORD_ANGSTROM:
          multiplier = 1.0;
          break;
          case QE_COORD_BOHR:
          model->fractional = FALSE;
          multiplier = BOHR_TO_ANGS;
          break;
          case QE_COORD_FRAC:
          break;
        }
      for (i=0; i< 3; i++)
        {
        fgetqeline(fp, line);
        buff = tokenize(line, &num_tokens);

        if (num_tokens == 3)
          {
#if DEBUG_READ_QE_OUT
          printf("cell params %s %s %s\n", buff[0], buff[1], buff[2]);
#endif
          model->latmat[0+i] = str_to_float(*(buff+0))*multiplier;
          model->latmat[3+i] = str_to_float(*(buff+1))*multiplier;
          model->latmat[6+i] = str_to_float(*(buff+2))*multiplier;

          g_strfreev(buff);
          }
        else
          {
          g_strfreev(buff);
          break;
          }
        }
      }
    }
  
  /* sort out lattice vectors */
  model->construct_pbc = TRUE;
  celldm[0] *= BOHR_TO_ANGS;
  switch (ibrav)
  {
    case 0:
    break;
    case 1:
    model->latmat[0] = model->latmat[4] = model->latmat[8] = celldm[0];
    break;
    case 4:
    model->latmat[0] = celldm[0];
    model->latmat[1] = celldm[0]*-0.5;
    model->latmat[4] = celldm[0]*0.5*sqrt(3);
    model->latmat[8] = celldm[0]*celldm[2];
    break;
    default:
    text = g_strdup_printf("ibrav = %d not yet implemented\n", ibrav);
    gui_text_show(ERROR, text);
    g_free(text);
    return(2);
  }
  
  strcpy(model->filename, filename);
  g_free(model->basename);
  model->basename = parse_strip(filename);
  
  model_prep(model);
  
  return(0);
}

