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
#include <stdlib.h>
#include <gmodule.h>

#include "gdis.h"
#include "file.h"
#include "parse.h"
#include "task.h"
#include "opengl.h"
#include "numeric.h"
#include "interface.h"

/* main data structures */
extern struct sysenv_pak sysenv;

/* data */
struct module_pak
{
/* location */
gchar *path;
/* top level description */
gchar *label;
/* make resident TRUE/FALSE */
gint resident;
/* list of gchar's */
GSList *symbols;
};

/*****************************/
/* free all modules and data */
/*****************************/
void module_free(void)
{
GSList *list;
struct module_pak *module;

for (list=sysenv.module_list ; list ; list=g_slist_next(list))
  {
  module = list->data;

  g_free(module->path);
  g_free(module->label);
  free_slist(module->symbols);
  g_free(module);
  }
g_slist_free(sysenv.module_list);
sysenv.module_list = NULL;
}

/****************************/
/* module config primitives */
/****************************/
void module_label_set(gchar *label, gpointer data)
{
struct module_pak *module = data;

module->label = g_strdup(label);
}

void module_symbol_add(gchar *symbol, gpointer data)
{
struct module_pak *module = data;

module->symbols = g_slist_prepend(module->symbols, g_strdup(symbol));
}

void module_resident_set(gint resident, gpointer data)
{
struct module_pak *module = data;

module->resident = resident;
}

/********************************/
/* module extraction primitives */
/********************************/
gchar *module_label(gpointer gp_module)
{
struct module_pak *module = gp_module;

return(module->label);
}

GSList *module_symbols(gpointer gp_module)
{
struct module_pak *module = gp_module;

return(module->symbols);
}

/*****************************************/
/* test if module is a valid gdis module */
/*****************************************/
#define DEBUG_MODULE_CONFIG 1
gpointer module_config(const gchar *path, const gchar *name)
{
gchar *basename, *fullname, *symbol;
GModule *gmodule;
gpointer ptr=NULL;
gint (*function)(gpointer);
struct module_pak *module=NULL;

/* checks */
if (!g_module_supported())
  {
  gui_text_show(WARNING, "Dynamic loading of modules is unsupported.\n");
  return(NULL);
  }

/* attempt to load module */
fullname = g_module_build_path(path, name);

#if DEBUG_MODULE_CONFIG
printf("Attempting to load: %s\n", fullname);
#endif

gmodule = g_module_open(fullname, 0);
if (gmodule)
  {
  basename = parse_strip(name);
  symbol = g_strdup_printf("module_%s_config", basename);
  if (g_module_symbol(gmodule, symbol, &ptr))
    {
/* attempt to configure module */
    function = ptr;
    module = g_malloc(sizeof(struct module_pak));
    module->symbols = NULL;
    if (function(module))
      {
#if DEBUG_MODULE_CONFIG
printf("ERROR: module configure function returned error.\n");
#endif
      g_free(module);
      module = NULL;
      }
    else
      {
/* sucessful configure */
      module->path = g_strdup(fullname);
      sysenv.module_list = g_slist_prepend(sysenv.module_list, module);
/* prevent module from being unloaded */
      if (module->resident)
        g_module_make_resident(gmodule);
      }
    }
  else
    {
#if DEBUG_MODULE_CONFIG
printf("ERROR: failed to locate config symbol: [%s]\n", symbol);
#endif
    ptr = NULL;
    }
  g_free(basename);
  g_module_close(gmodule);
  }
#if DEBUG_MODULE_CONFIG
else
  printf("%s\n", g_module_error());
#endif

/* cleanup */
g_free(fullname);
return(ptr);
}

/**************************/
/* scan for valid modules */
/**************************/
#define DEBUG_MODULE_SETUP 1
void module_setup(void)
{
gchar *name, *ext;
GSList *dir, *list;

/* search gdis directory for valid plug-ins */
dir = file_dir_list(sysenv.gdis_path, FALSE);
for (list=dir ; list ; list=g_slist_next(list))
  {
  name = list->data;

  ext = file_extension_get(name);

  if (ext)
    {
    if (strlen(ext) == 2)
      {
      if (g_ascii_strncasecmp(ext, "so", 2) == 0)
        {
        if (module_config(sysenv.gdis_path, name))
          {
#if DEBUG_MODULE_SETUP
printf("module: %s is valid.\n", name);
#endif
          }
        else
          {
#if DEBUG_MODULE_SETUP
printf("module: %s is invalid.\n", name);
#endif
          }
        }
      }
    }
  }
free_slist(dir);
}

/*********************/
/* module invocation */
/*********************/
gint module_invoke(gchar *name, gpointer arg, gpointer ptr_module)
{
gint value = 0;
GModule *gmodule;
gpointer symbol;
void (*function)(gpointer);
struct module_pak *module;

/* attempt to open module */
module = ptr_module;
gmodule = g_module_open(module->path, 0);
if (gmodule)
  {
  if (!g_module_symbol(gmodule, name, &symbol))
    value = 2;
  }
else
  value = 1;

/* invoke module function */
if (!value)
  {
  function = symbol;
  function(arg);
  }  

/* done */
g_module_close(gmodule);
return(value);
}

