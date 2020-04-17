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

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "graph.h"
#include "surface.h"
#include "molsurf.h"
#include "numeric.h"
#include "parse.h"
#include "project.h"
#include "file.h"
#include "scan.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/***************************/
/* destroy a configuration */
/***************************/
void project_config_free(gpointer data)
{
}

/**********************/
/* destroy a projects */
/**********************/
void project_free(gpointer data)
{
struct project_pak *project = data;

g_assert(project != NULL);

g_free(project->label);
g_free(project->path);
g_slist_free(project->models);
g_slist_free(project->configs);

/* TODO - free (project specific) hash table data */
g_hash_table_destroy(project->data);

g_free(project);
}

/************************/
/* create a new project */
/************************/
gpointer project_new(const gchar *label)
{
struct project_pak *project;

project = g_malloc(sizeof(struct project_pak));
sysenv.projects = g_slist_append(sysenv.projects, project);

project->path = NULL;
if (label)
  project->label = g_strdup(label);
else
  project->label = g_strdup("project x");

project->models = NULL;
project->configs = NULL;
project->data = g_hash_table_new(g_str_hash, g_str_equal);

return(project);
}

/****************************/
/* add a model to a project */
/****************************/
void project_model_add(struct model_pak *model, gpointer data)
{
struct project_pak *project = data;

project->models = g_slist_append(project->models, model);
}

/***********************************/
/* create a new model in a project */
/***********************************/
gpointer project_model_new(gpointer data)
{
struct project_pak *project = data;
struct model_pak *model;

/* NB: too many things look at model_pak -> keep public */
model = g_malloc(sizeof(struct model_pak));
project->models = g_slist_append(project->models, model);

model_init(model);
model->project = project;

return(model);
}

/***************************/
/* project data primitives */
/***************************/
void project_data_set(gchar *key, gpointer data, gpointer p)
{
struct project_pak *project = p;

g_assert(project != NULL);

g_hash_table_insert(project->data, key, data);
}

gpointer project_data_get(gchar *key, gpointer p)
{
struct project_pak *project = p;

g_assert(project != NULL);

return(g_hash_table_lookup(project->data, key));
}

/***********************************************/
/* create a new computional code configuration */
/***********************************************/
gpointer project_config_get(gint id, gpointer data)
{
gchar *label;
GSList *list;
struct config_pak *config;
struct project_pak *project = data;

g_assert(project != NULL);

/* search for pre-existing config */
for (list=project->configs ; list ; list=g_slist_next(list))
  {
  config = list->data;
  if (config->id == id)
    return(config);
  }

/* create new config */
switch (id)
  {
  case ABINIT:
    label = g_strdup("abinit");
    break;
  case CASTEP:
    label = g_strdup("castep");
    break;
  case GAMESS:
    label = g_strdup("gamess");
    break;
  case GAUSS:
    label = g_strdup("gaussian");
    break;
  case GULP:
    label = g_strdup("gulp");
    break;

  default:
    gui_text_show(ERROR, "Unknown computational code.\n");
    return(NULL);
  }

config = g_malloc(sizeof(struct config_pak));
/* TODO - add using a sort */
project->configs = g_slist_prepend(project->configs, config);
config->id = id;
config->label = label;
config->data = NULL;

return(config);
}

/************************/
/* extract project path */
/************************/
gchar *project_path(gpointer data)
{
struct project_pak *project = data;

g_assert(project != NULL);

return(project->path);
}

/****************************************************/
/* extract a pointer to the nth model of a projects */
/****************************************************/
gpointer project_model_get(gint n, gpointer data)
{
struct project_pak *project = data;

g_assert(project != NULL);

return(g_slist_nth_data(project->models, n));
}

/*************************/
/* process solvent data */
/*************************/
void project_solvent_plot(GArray *energy, struct model_pak *model)
{
gint i, j, bins;
gdouble value, min, max, de, *plot;
gpointer graph;

g_assert(model != NULL);
g_assert(energy != NULL);

if (!energy->len)
  return;

min = max = g_array_index(energy, gdouble, 0);
for (i=0 ; i<energy->len ; i++)
  {
  value = g_array_index(energy, gdouble, i);
  
  if (value < min)
    min = value;
  if (value > max)
    max = value;
  }

de = max - min;

bins = 50;
plot = g_malloc(bins * sizeof(gdouble));
for (i=bins ; i-- ; )
  plot[i] = 0.0;

for (i=0 ; i<energy->len ; i++)
  {
  value = g_array_index(energy, gdouble, i);
  value -= min;
  value /= de;
  value *= bins;
  j = (gint) (value - 0.5);

  j = CLAMP(j, 0, bins-1);

  plot[j]++;
  }

graph = graph_new("Solvent", model);
graph_add_data(bins, plot, min, max, graph);

g_free(plot);
}

/**********************************/
/* read in a project control file */
/**********************************/
#define DEBUG_PROJECT_READ 0
gint project_read(gchar *filename, struct model_pak *model)
{
gint a, b, i, j, ii, ij, ix, jx, ia, ib;
gint set, fail, type, size, num_tokens;
gint grid[2];
gchar *basename, *fullname, **buff;
gdouble e, min, max, x[2], cell[2], *surface=NULL;
gpointer scan, project;
GArray *energy;
GSList *list, *list_bad=NULL, *list_error=NULL, *list_nomin=NULL;
struct model_pak temp_model;

g_assert(model != NULL);
g_assert(filename != NULL);

/* init */
fullname = g_build_filename(sysenv.cwd, filename, NULL);
model_init(&temp_model);
temp_model.grafted = TRUE;

grid[0]=0;
grid[1]=0;

/* CURRENT */
project = project_new("solvent");
project_model_add(model, project);
model->project = project;

#if DEBUG_PROJECT_READ
printf("Reading project file: %s\n", fullname);
#endif

scan = scan_new(fullname);
g_free(fullname);

energy = g_array_new(FALSE, FALSE, sizeof(gdouble));

a = b = ia = ib = size = set = type = 0;
while (!scan_complete(scan))
  {
  buff = scan_get_tokens(scan, &num_tokens);

/* keyword checks */
  if (buff)
    {
    if (g_ascii_strncasecmp(*buff, "%set", 4) == 0)
      {
      switch (type)
        {
        case 1:
          if (num_tokens == 5)
            {
/* grid points */
            grid[0] = str_to_float(*(buff+1));
            grid[1] = str_to_float(*(buff+2));
/* sampling fraction */
            cell[0] = str_to_float(*(buff+3));
            cell[1] = str_to_float(*(buff+4));
/* total number of grid points on surface (includes symmetry related images) */
            a = nearest_int(grid[0]/cell[0]);
            b = nearest_int(grid[1]/cell[1]);
/* number of images required to fill entire cell */
            ia = nearest_int(1.0/cell[0]);
            ib = nearest_int(1.0/cell[1]);
/* allocation */
            size = a*b;
	    if(surface!=NULL) g_free(surface);/*FIX a02228*/
            surface = g_malloc(size*sizeof(gdouble));
            for (i=size ; i-- ; )
              surface[i] = 0.0;
            set = 1;
            }
          else
            {
printf("Bad solvent set.\n");
g_assert_not_reached();
            }
          break;

        default:
          set = 1;
        }
      g_strfreev(buff);
      continue;
      }

    if (g_ascii_strncasecmp(*buff, "%title", 6) == 0)
      {
/* FIXME - properly defined project types */
      if (num_tokens > 1)
        {
        if (g_ascii_strncasecmp(*(buff+1), "solvent", 7) == 0)
          type = 1;
        }
      g_strfreev(buff);
      continue;
      }
    }
  else
    {
    //set=0;/*FIX 11b326*/
    break;
    }

/* processing */
  if (set)
    {
/*
printf("set dimension: %d x %d, scaling: %f x %f\n", max[0], max[1], scale[0], scale[1]);
*/

    basename = parse_strip(*buff);
    fullname =  g_strdup_printf("%s%s%s.got", sysenv.cwd, DIR_SEP, basename);
    fail = 0;

#if DEBUG_PROJECT_READ
printf("reading: %s      ", fullname);
#endif

/* bad read check */
    temp_model.error_file_read = g_string_assign(temp_model.error_file_read, "");
    if (read_gulp_output(fullname, &temp_model))
      {
      list_bad = g_slist_prepend(list_bad, g_strdup_printf("%s.got", basename));
      fail++;
      }

/* check for any GULP errors */
    if ((temp_model.error_file_read)->len)
      {
      list_error = g_slist_prepend(list_error, g_strdup_printf("%s.got", basename));
      fail++;
      }

/* ignore anything that didn't converge */
    if (temp_model.gulp.gnorm > 0.1)
      {
      list_nomin = g_slist_prepend(list_nomin, g_strdup_printf("%s.got", basename));
      fail++;
      }

/* we have a valid energy that can be processed */
    if (!fail)
     {
#if DEBUG_PROJECT_READ
printf("[%f][%f]\n", temp_model.gulp.energy, temp_model.gulp.gnorm);
#endif

     g_array_append_val(energy, temp_model.gulp.energy);

/* project dependent parsing */
     if (surface)
       {
       switch (type)
         {
         case 1:
/* grid point -> array index */
           x[0] = str_to_float(*(buff+1));
           x[1] = str_to_float(*(buff+2));
           x[0] *= a;
           x[1] *= b;
           i = nearest_int(x[0]);
           j = nearest_int(x[1]);
/* images (if any) */
           for (ij=0 ; ij<ib ; ij++)
             {
             for (ii=0 ; ii<ia ; ii++)
               {
               ix = i + ii*grid[0];
               jx = j + ij*grid[1];
               if (temp_model.gulp.energy < *(surface+jx*a+ix))
                 *(surface+jx*a+ix) = temp_model.gulp.energy;
               }
             }
            break;
          }
        }
      }
    g_free(fullname);
    }
  g_strfreev(buff);
  }
scan_free(scan);

/* output summary for project import */
printf("=======================================\n");
if (list_bad)
  {
  list_bad = g_slist_reverse(list_bad);
  printf("The following files were bad or missing:\n");
  for (list=list_bad ; list ; list=g_slist_next(list))
    {
    printf("%s\n", (gchar *) list->data);
    }
  printf("=======================================\n");
  }
if (list_error)
  {
  list_error = g_slist_reverse(list_error);
  printf("The following files had GULP errors:\n");
  for (list=list_error ; list ; list=g_slist_next(list))
    {
    printf("%s\n", (gchar *) list->data);
    }
  printf("=======================================\n");
  }
if (list_nomin)
  {
  list_nomin = g_slist_reverse(list_nomin);
  printf("The following files failed to minimize:\n");
  for (list=list_nomin ; list ; list=g_slist_next(list))
    {
    printf("%s\n", (gchar *) list->data);
    }
  printf("=======================================\n");
  }

/* process the data set(s) */
switch (type)
  {
  case 1:
#if DEBUG_PROJECT_READ
printf("Processing solvent data...\n");
#endif

    min = G_MAXDOUBLE;
    max = -G_MAXDOUBLE;
    for (i=size ; i-- ; )
      {
      e = *(surface+i);

      if (e < min)
        min = e;
      if (e > max)
        max = e;
      }

/* setup colour scale */
/* TODO - make the names generic */
    model->epot_min = min;
    model->epot_max = max;
    model->epot_div = 11;

/* plot molsurf and colour with dock energy */
    project_data_set("dock:a", GINT_TO_POINTER(a), project);
    project_data_set("dock:b", GINT_TO_POINTER(b), project);

    project_data_set("dock:min", g_strdup_printf("%f", min), project);
    project_data_set("dock:max", g_strdup_printf("%f", max), project);

    project_data_set("dock:energy", surface, project);
    ms_cube(0.8, MS_MOLECULAR, MS_SOLVENT, model);
    coords_init(CENT_COORDS, model);

/* graph the solvent interaction */
    project_solvent_plot(energy, model);
    break;

  default:
    printf("Unknown project type.\n");
    break;
  }

/* cleanup */
model_free(&temp_model);
g_array_free(energy, TRUE);
tree_model_refresh(model);
free_slist(list_bad);
free_slist(list_error);
free_slist(list_nomin);

return(0);
}
