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
#include <strings.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "file.h"
#include "model.h"
#include "parse.h"
#include "render.h"
#include "spatial.h"
#include "matrix.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* triggers */
enum 
{
XML_INACTIVE,
XML_SYSTEM, XML_SYMMETRY, XML_CRYSTAL, XML_CELL, XML_CORE, XML_SHELL,
XML_CELL_A, XML_CELL_B, XML_CELL_C,
XML_CELL_ALPHA, XML_CELL_BETA, XML_CELL_GAMMA,
XML_WAYPOINT, XML_SPATIAL, XML_VERTEX, XML_PICTURE
};

gint xml_context = XML_INACTIVE;
struct model_pak *xml_model=NULL;
struct core_pak *xml_core=NULL;
struct shel_pak *xml_shell=NULL;
GHashTable *xml_system_table;


/* TEMP - parse diffax sfc data */
void parse_sfc(FILE *fp)
{
gint i, n;
gchar *line, *type, **buff;
gdouble *sfc;
GSList *list;

while ((line = file_read_line(fp)))
  {
  type = g_strstrip(g_strndup(line, 6));

  if (elem_symbol_test(type))
    {
/* tokenize everything after the atom type */
    buff = tokenize(&line[6], &n);

printf("[%s] ", type);
    list = NULL;
    for (i=0 ; i<n ; i++)
      {
      sfc = g_malloc(sizeof(gdouble));
      *sfc = str_to_float(*(buff+i));

printf("[%f] ", *sfc);
      list = g_slist_prepend(list, sfc);
      }
    list = g_slist_reverse(list);
printf("\n");

    g_hash_table_insert(sysenv.sfc_table, type, list);

    g_strfreev(buff);
    }
  g_free(line);
  }
}


/************************************/
/* find or allocate a model pointer */
/************************************/
void xml_parse_system(const gchar **names, const gchar **values)
{
gint i;
gchar *key;

/* get the system id value */
i=0;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "id\0", 3) == 0)
    {
/* first time exception - model already allocated */
    if (xml_model)
      {
printf("initial system: %s\n", *(values+i));
      key = g_strdup(*(values+i));
      g_hash_table_insert(xml_system_table, key, xml_model);
      }
    else
      {
/* check id against existing models */
      xml_model = g_hash_table_lookup(xml_system_table, *(values+i)); 
      if (!xml_model)
        {
printf("adding system: %s\n", *(values+i));
/* id not found - allocate new model & add to hash table */
        xml_model = model_new();
        key = g_strdup(*(values+i));
        g_hash_table_insert(xml_system_table, key, xml_model);
        }
      else
printf("existing system: %s\n", *(values+i));
      }
    }
  i++;
  }
/* got all we want */
xml_context = XML_INACTIVE;
}

/***************************************/
/* process symmetry element attributes */
/***************************************/
void xml_parse_symmetry(const gchar **names, const gchar **values)
{
gint i;

g_assert(xml_model != NULL);

/* process attributes */
i=0;
while (*(names+i))
  {
/* only init the cell parsing trigger if it's the full cell's data */
  if (g_ascii_strncasecmp(*(names+i), "spacegroup", 10) == 0)
    {
    xml_model->sginfo.spacename = g_strdup(*(values+i));
    }
  i++;
  }
}

/********************************************/
/* start cell reading if we get a full cell */
/********************************************/
void xml_parse_fullcell(const gchar **names, const gchar **values)
{
gint i;

/* process attributes */
i=0;
xml_context = XML_INACTIVE;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "dictRef\0", 8) == 0)
    {
    if (g_ascii_strncasecmp(*(values+i), "gulp:fullcell\0", 16) == 0)
      {
      xml_context = XML_CELL;
      xml_model->periodic = 3;
      }
    if (g_ascii_strncasecmp(*(values+i), "gulp:surface\0", 13) == 0)
      {
      xml_context = XML_CELL;
      xml_model->periodic = 2;
      }
    if (g_ascii_strncasecmp(*(values+i), "gulp:polymer\0", 13) == 0)
      {
      xml_context = XML_CELL;
      xml_model->periodic = 1;
      }
    }
  i++;
  }
}

/***************************/
/* process cell attributes */
/***************************/
void xml_parse_cell(const gchar **names, const gchar **values)
{
gint i;

g_assert(xml_model != NULL);

/* process attributes (if any) */
i=0;
while (*(names+i))
  {
/* TODO - CHECK that the units are angs and degrees */
  if (g_ascii_strncasecmp(*(names+i), "dictRef\0", 8) == 0)
    {
/* TODO - split (remove namespace) and check a/b/etc. only? */
/* cell value require the text callback */
    if (g_ascii_strncasecmp(*(values+i), "cml:a\0", 6) == 0)
      xml_context = XML_CELL_A;
    if (g_ascii_strncasecmp(*(values+i), "cml:b\0", 6) == 0)
      xml_context = XML_CELL_B;
    if (g_ascii_strncasecmp(*(values+i), "cml:c\0", 6) == 0)
      xml_context = XML_CELL_C;
    if (g_ascii_strncasecmp(*(values+i), "cml:alpha\0", 10) == 0)
      xml_context = XML_CELL_ALPHA;
    if (g_ascii_strncasecmp(*(values+i), "cml:beta\0", 9) == 0)
      xml_context = XML_CELL_BETA;
    if (g_ascii_strncasecmp(*(values+i), "cml:gamma\0", 10) == 0)
      xml_context = XML_CELL_GAMMA;
    }
  i++;
  }
}

/*******************/
/* parse atom info */
/*******************/
void xml_parse_atom(const gchar **names, const gchar **values)
{
gint i, n;
gchar **buff;

g_assert(xml_model != NULL);

/* init structure */
if (!xml_core)
  {
  xml_core = new_core("X", xml_model);
  xml_model->cores = g_slist_prepend(xml_model->cores, xml_core);
  }

/* process attributes (if any) */
i=0;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "elementType\0", 12) == 0)
    {
/* FIXME - potential overflow bug here, since label is an array */
    g_free(xml_core->atom_label);
    xml_core->atom_label = g_strdup(*(values+i));     
    core_init(xml_core->atom_label, xml_core, xml_model);
    }
  if (g_ascii_strncasecmp(*(names+i), "xyzf", 4) == 0)
    {
    buff = tokenize(*(values+i), &n);
    if (n > 2)
      {
      xml_core->x[0] = str_to_float(*(buff+0));
      xml_core->x[1] = str_to_float(*(buff+1));
      xml_core->x[2] = str_to_float(*(buff+2));
      xml_model->fractional = TRUE;
      }
    }
  if (g_ascii_strncasecmp(*(names+i), "xf", 2) == 0)
    {
    xml_core->x[0] = str_to_float(*(values+i));
/* FIXME - inefficient to repeat this for all atoms */
/* but may be necessary if every atom is allowed to be cart or fract */
    xml_model->fractional = TRUE;
    }
  if (g_ascii_strncasecmp(*(names+i), "x3", 2) == 0)
    {
    xml_core->x[0] = str_to_float(*(values+i));
/* FIXME - inefficient to repeat this for all atoms */
/* but may be necessary if every atom is allowed to be cart or fract */
    xml_model->fractional = FALSE;
    }

  if (g_ascii_strncasecmp(*(names+i), "yf", 2) == 0 ||
      g_ascii_strncasecmp(*(names+i), "y3", 2) == 0) 
    {
    xml_core->x[1] = str_to_float(*(values+i));
    }
  if (g_ascii_strncasecmp(*(names+i), "zf", 2) == 0 ||
      g_ascii_strncasecmp(*(names+i), "z3", 2) == 0) 
    {
    xml_core->x[2] = str_to_float(*(values+i));
    }
  i++;
  }
}

/********************/
/* parse shell info */
/********************/
void xml_parse_shell(const gchar **names, const gchar **values)
{
gint i, n;
gchar **buff;

g_assert(xml_model != NULL);
g_assert(xml_core != NULL);

/* init structure */
if (!xml_shell)
  {
  xml_shell = new_shell(xml_core->atom_label, xml_model);
  xml_model->shels = g_slist_prepend(xml_model->shels, xml_shell);
  }

/* process attributes (if any) */
i=0;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "xyzf", 4) == 0)
    {
    buff = tokenize(*(values+i), &n);
    if (n > 2)
      {
      xml_core->x[0] = str_to_float(*(buff+0));
      xml_core->x[1] = str_to_float(*(buff+1));
      xml_core->x[2] = str_to_float(*(buff+2));
      xml_model->fractional = TRUE;
      }
    }
  if (g_ascii_strncasecmp(*(names+i), "xf", 2) == 0)
    {
    xml_core->x[0] = str_to_float(*(values+i));
/* FIXME - inefficient to repeat this for all atoms */
/* but may be necessary if every atom is allowed to be cart or fract */
    xml_model->fractional = TRUE;
    }
  if (g_ascii_strncasecmp(*(names+i), "x3", 2) == 0)
    {
    xml_core->x[0] = str_to_float(*(values+i));
/* FIXME - inefficient to repeat this for all atoms */
/* but may be necessary if every atom is allowed to be cart or fract */
    xml_model->fractional = FALSE;
    }

  if (g_ascii_strncasecmp(*(names+i), "yf", 2) == 0 ||
      g_ascii_strncasecmp(*(names+i), "y3", 2) == 0) 
    {
    xml_core->x[1] = str_to_float(*(values+i));
    }
  if (g_ascii_strncasecmp(*(names+i), "zf", 2) == 0 ||
      g_ascii_strncasecmp(*(names+i), "z3", 2) == 0) 
    {
    xml_core->x[2] = str_to_float(*(values+i));
    }
  i++;
  }
}

/*******************/
/* spatial parsing */
/*******************/
struct spatial_pak *xml_spatial;

/************************/
/* create a new spatial */
/************************/
void xml_parse_spatial(const gchar **names, const gchar **values)
{
gint i, periodic=-1;

g_assert(xml_model != NULL);

i=0;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "periodic", 8) == 0)
    periodic = str_to_float(*(values+i));
/* TODO - more details */
  i++;
  }

xml_spatial = spatial_new(NULL, SPATIAL_GENERIC, 0, periodic, xml_model);
}

/***************************/
/* create camera waypoints */
/***************************/
void xml_parse_waypoint(const gchar **names, const gchar **values)
{
gint i;
struct camera_pak *camera;

g_assert(xml_model != NULL);

camera = camera_new();

i=0;
while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "mode", 4) == 0)
    camera->mode = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "perspective", 11) == 0)
    camera->perspective = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "fov", 3) == 0)
    camera->fov = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "zoom", 4) == 0)
    camera->zoom = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "x0", 2) == 0)
    camera->x[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "x1", 2) == 0)
    camera->x[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "x2", 2) == 0)
    camera->x[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "o0", 2) == 0)
    camera->o[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "o1", 2) == 0)
    camera->o[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "o2", 2) == 0)
    camera->o[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "v0", 2) == 0)
    camera->v[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "v1", 2) == 0)
    camera->v[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "v2", 2) == 0)
    camera->v[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "e0", 2) == 0)
    camera->e[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "e1", 2) == 0)
    camera->e[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "e2", 2) == 0)
    camera->e[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "q0", 2) == 0)
    camera->q[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "q1", 2) == 0)
    camera->q[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "q2", 2) == 0)
    camera->q[2] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "q3", 2) == 0)
    camera->q[3] = str_to_float(*(values+i));

  i++;
  }
xml_model->waypoint_list = g_slist_prepend(xml_model->waypoint_list, camera);
}

/********************/
/* populate spatial */
/********************/
void xml_parse_vertex(const gchar **names, const gchar **values)
{
gint i=0;
gdouble x[3], n[3], c[3];

g_assert(xml_spatial != NULL);

VEC3SET(x, 0.0, 0.0, 0.0);
VEC3SET(n, 0.0, 0.0, 0.0);
VEC3SET(c, 0.0, 0.0, 0.0);

while (*(names+i))
  {
  if (g_ascii_strncasecmp(*(names+i), "red", 3) == 0)
    c[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "green", 5) == 0)
    c[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "blue", 4) == 0)
    c[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "x3", 2) == 0)
    x[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "y3", 2) == 0)
    x[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "z3", 2) == 0)
    x[2] = str_to_float(*(values+i));

  if (g_ascii_strncasecmp(*(names+i), "nx", 2) == 0)
    n[0] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "ny", 2) == 0)
    n[1] = str_to_float(*(values+i));
  if (g_ascii_strncasecmp(*(names+i), "nz", 2) == 0)
    n[2] = str_to_float(*(values+i));

  i++;
  }
spatial_vnorm_add(x, n, c, xml_spatial);
}

/****************/
/* process text */
/****************/
void xml_parse_text(GMarkupParseContext *context,
                    const gchar *text,
                    gsize text_len,  
                    gpointer user_data,
                    GError **error)
{
switch (xml_context)
  {
/* specific cell contexts */
  case XML_CELL_A:
    xml_model->pbc[0] = str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;
  case XML_CELL_B:
    xml_model->pbc[1] = str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;
  case XML_CELL_C:
    xml_model->pbc[2] = str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;
  case XML_CELL_ALPHA:
    xml_model->pbc[3] = D2R * str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;
  case XML_CELL_BETA:
    xml_model->pbc[4] = D2R * str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;
  case XML_CELL_GAMMA:
    xml_model->pbc[5] = D2R * str_to_float(text);
/* revert to general cell value context */
    xml_context = XML_CELL;
    break;

  case XML_PICTURE:
    xml_model->picture_list = g_slist_append(xml_model->picture_list,
                                         g_strstrip(g_strdup(text)));
    break;
  }
}

/**************************/
/* element start callback */
/**************************/
void xml_start_element(GMarkupParseContext *context,
                       const gchar *element_name,
                       const gchar **attribute_names,
                       const gchar **attribute_values,
                       gpointer current_model,
                       GError **error)
{
struct model_pak *model;

model = current_model;

/* context switching control */
if (g_ascii_strncasecmp(element_name, "system\0", 7) == 0)
  xml_context = XML_SYSTEM;
if (g_ascii_strncasecmp(element_name, "crystal\0", 8) == 0)
  xml_context = XML_CRYSTAL;
if (g_ascii_strncasecmp(element_name, "symmetry\0", 9) == 0)
  xml_context = XML_SYMMETRY;
if (g_ascii_strncasecmp(element_name, "atom\0", 5) == 0)
  xml_context = XML_CORE;
if (g_ascii_strncasecmp(element_name, "particle\0", 9) == 0)
  xml_context = XML_SHELL;
if (g_ascii_strncasecmp(element_name, "spatial\0", 8) == 0)
  xml_context = XML_SPATIAL;
if (g_ascii_strncasecmp(element_name, "vertex\0", 7) == 0)
  xml_context = XML_VERTEX;
if (g_ascii_strncasecmp(element_name, "picture\0", 8) == 0)
  xml_context = XML_PICTURE;
if (g_ascii_strncasecmp(element_name, "waypoint\0", 9) == 0)
  xml_context = XML_WAYPOINT;

/* context attribute parsing */
switch (xml_context)
  {
  case XML_SYSTEM:
    xml_parse_system(attribute_names, attribute_values);
    return;

  case XML_CRYSTAL:
    xml_parse_fullcell(attribute_names, attribute_values);
    return;

  case XML_CELL:
    xml_parse_cell(attribute_names, attribute_values);
    return;

  case XML_CORE:
    xml_parse_atom(attribute_names, attribute_values);
    return;

  case XML_SHELL:
    xml_parse_shell(attribute_names, attribute_values);
    return;

  case XML_SYMMETRY:
    xml_parse_symmetry(attribute_names, attribute_values);
    return;

  case XML_SPATIAL:
    xml_parse_spatial(attribute_names, attribute_values);
    return;

  case XML_VERTEX:
    xml_parse_vertex(attribute_names, attribute_values);
    return;

  case XML_WAYPOINT:
    xml_parse_waypoint(attribute_names, attribute_values);
    return;
  }
}

/************************/
/* element end callback */
/************************/
void xml_end_element(GMarkupParseContext *context,
                     const gchar *element_name,
                     gpointer current_model,
                     GError **error)
{
/* reset parsing context and related data */
if (g_ascii_strncasecmp(element_name, "system\0", 7) == 0)
  {
/* end of a model - prep? */
  if (xml_model->waypoint_list)
    xml_model->waypoint_list = g_slist_reverse(xml_model->waypoint_list);

  xml_model = NULL;
  xml_context = XML_INACTIVE;
  }
if (g_ascii_strncasecmp(element_name, "crystal\0", 8) == 0)
  {
  xml_context = XML_INACTIVE;
  }
if (g_ascii_strncasecmp(element_name, "atom\0", 5) == 0)
  {
  xml_core = NULL;
  xml_context = XML_INACTIVE;
  }
if (g_ascii_strncasecmp(element_name, "particle\0", 5) == 0)
  {
  xml_shell = NULL;
  xml_context = XML_CORE;
  }
if (g_ascii_strncasecmp(element_name, "spatial\0", 8) == 0)
  {
  xml_context = XML_INACTIVE;
  xml_spatial->list = g_slist_reverse(xml_spatial->list);
  xml_spatial = NULL;
  }
if (g_ascii_strncasecmp(element_name, "vertex\0", 7) == 0)
  xml_context = XML_SPATIAL;
if (g_ascii_strncasecmp(element_name, "picture\0", 8) == 0)
  xml_context = XML_INACTIVE;
if (g_ascii_strncasecmp(element_name, "waypoint\0", 9) == 0)
  xml_context = XML_INACTIVE;
}

/******************************/
/* read a single model's data */
/******************************/
gint read_xml_frame(FILE *fp, struct model_pak *model)
{
gchar *line;
GMarkupParser xml_parser;
GMarkupParseContext *xml_context;
GError *error;

xml_system_table = g_hash_table_new_full(&g_str_hash, &hash_strcmp, &g_free, NULL);

/* TODO - think more about this... */
xml_model = model;

/* setup context parse (ie callbacks) */
xml_parser.start_element = &xml_start_element;
xml_parser.end_element = &xml_end_element;
xml_parser.text = &xml_parse_text;
xml_parser.passthrough = NULL;
xml_parser.error = NULL;
xml_context = g_markup_parse_context_new(&xml_parser, 0, model, NULL);

/* read in blocks (lines) of text */
line = file_read_line(fp);
while (line)
  {
/* parse the line */
  if (!g_markup_parse_context_parse(xml_context, line, strlen(line), &error))
    printf("read_xml() : parsing error.\n");

  g_free(line);
  line = file_read_line(fp);
  }

/* cleanup */
if (!g_markup_parse_context_end_parse(xml_context, &error))
  printf("read_xml() : errors occurred reading file.\n");

g_markup_parse_context_free(xml_context);
g_hash_table_destroy(xml_system_table);

return(0);
}

/*****************************/
/* the main reading function */
/*****************************/
#define DEBUG_READ_XML 0
gint read_xml(gchar *filename, struct model_pak *model)
{
gint /*i, n,*/ num_models=0;
FILE *fp;

fp = fopen(filename, "rt");
if (!fp)
  return(1);

/* TODO - store current file position in frame list */
/* TODO - alloc new model pointers for the read */
read_xml_frame(fp, model);
num_models++; 

fclose(fp);

#if DEBUG_READ_XML
printf("read_xml() : found %d models.\n", num_models);
#endif

model_prep(model);

return(0);
}

/*****************************/
/* the main writing function */
/*****************************/
gint write_xml(gchar *filename, struct model_pak *model)
{
gdouble vec[3];
GSList *list, *list2;
struct core_pak *core;
struct vec_pak *v;
struct spatial_pak *spatial;
FILE *fp;

fp = fopen(filename, "wt");
if (!fp)
  return(1);

fprintf(fp, "<system>\n");

/* periodicity */
if (model->periodic)
  {
  switch (model->periodic)
    {
    case 1:
      fprintf(fp, "  <crystal dictRef=\"gulp:polymer\">\n");
      break;

    case 2:
      fprintf(fp, "  <crystal dictRef=\"gulp:surface\">\n");
      break;

    default:
      fprintf(fp, "  <crystal dictRef=\"gulp:fullcell\">\n");
      break;
    }
  fprintf(fp, "    <cell dictRef=\"cml:a\"> %f </cell>\n", model->pbc[0]);
  fprintf(fp, "    <cell dictRef=\"cml:b\"> %f </cell>\n", model->pbc[1]);
  fprintf(fp, "    <cell dictRef=\"cml:c\"> %f </cell>\n", model->pbc[2]);
  fprintf(fp, "    <cell dictRef=\"cml:alpha\"> %f </cell>\n",
                                           R2D*model->pbc[3]);
  fprintf(fp, "    <cell dictRef=\"cml:beta\"> %f </cell>\n",
                                          R2D*model->pbc[4]);
  fprintf(fp, "    <cell dictRef=\"cml:gamma\"> %f </cell>\n",
                                           R2D*model->pbc[5]);
  fprintf(fp, "  </crystal>\n");
  }

/* cores */
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  ARR3SET(vec, core->x);
  vecmat(model->latmat, vec);
/*
  fprintf(fp, "  <atom elementType=\"%s\" x3=\"%f\" y3=\"%f\" z3=\"%f\">",
          core->atom_label, vec[0], vec[1], vec[2]);
*/

  fprintf(fp, "  <atom elementType=\"%s\" x3=\"%f\" y3=\"%f\" z3=\"%f\">",
          elements[core->atom_code].symbol, vec[0], vec[1], vec[2]);


  fprintf(fp, "</atom>\n");
/* TODO - shells as particles */
  }

/* spatial data */
for (list=model->spatial ; list ; list=g_slist_next(list))
  {
  spatial = list->data;

  fprintf(fp, "<spatial Method=\"%d\" Periodic=\"%d\">\n",
                      spatial->method, spatial->periodic);

  for (list2=spatial->list ; list2 ; list2=g_slist_next(list2))
    {
    v = list2->data;

    fprintf(fp, "  <vertex red=\"%f\" green=\"%f\" blue=\"%f\"",
                 v->colour[0], v->colour[1], v->colour[2]);
 
    fprintf(fp, " x3=\"%f\" y3=\"%f\" z3=\"%f\"",
                 v->x[0], v->x[1], v->x[2]);

    fprintf(fp, " nx=\"%f\" ny=\"%f\" nz=\"%f\"></vertex>\n",
                 v->n[0], v->n[1], v->n[2]);
    }

  fprintf(fp, "</spatial>\n");
  }

/* camera waypoints */
for (list=model->waypoint_list ; list ; list=g_slist_next(list))
  {
  struct camera_pak *camera = list->data;

  fprintf(fp, "<waypoint mode=\"%d\" perspective=\"%d\" fov=\"%f\" zoom=\"%f\"\n",
              camera->mode, camera->perspective, camera->fov, camera->zoom);
  fprintf(fp, "x0=\"%f\" x1=\"%f\" x2=\"%f\"\n",camera->x[0],camera->x[1],camera->x[2]);
  fprintf(fp, "o0=\"%f\" o1=\"%f\" o2=\"%f\"\n",camera->o[0],camera->o[1],camera->o[2]);
  fprintf(fp, "v0=\"%f\" v1=\"%f\" v2=\"%f\"\n",camera->v[0],camera->v[1],camera->v[2]);
  fprintf(fp, "e0=\"%f\" e1=\"%f\" e2=\"%f\"\n",camera->e[0],camera->e[1],camera->e[2]);
  fprintf(fp, "q0=\"%f\" q1=\"%f\" q2=\"%f\" q3=\"%f\">", 
               camera->q[0], camera->q[1], camera->q[2], camera->q[3]);
  fprintf(fp, "</waypoint>\n"); 
  }

/* pictures */
for (list=model->picture_list ; list ; list=g_slist_next(list))
  fprintf(fp, "<picture> %s </picture>\n", (gchar *) list->data);

fprintf(fp, "</system>\n");

fclose(fp);
return(0);
}

