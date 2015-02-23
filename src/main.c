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

/* irix */
#define _BSD_SIGNALS 1

#include <signal.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#ifndef __WIN32
#include <sys/wait.h>
#endif

#include "gdis.h"
#include "command.h"
#include "file.h"
#include "parse.h"
#include "task.h"
#include "host.h"
#include "render.h"
#include "matrix.h"
#include "model.h"
#include "numeric.h"
#include "module.h"
#include "library.h"
#include "interface.h"
#include "gui_image.h"

/* main data structures */
struct sysenv_pak sysenv;
struct elem_pak elements[MAX_ELEMENTS];

/********************************/
/* sfc hash table value cleanup */
/********************************/
void free_sfc(gpointer sfc_list)
{
free_slist((GSList *) sfc_list);
}

/******************/
/* system cleanup */
/******************/
void sys_free(void)
{

file_free();

free_slist(sysenv.render.light_list);

g_hash_table_destroy(sysenv.sfc_table);
g_hash_table_destroy(sysenv.library);
g_hash_table_destroy(sysenv.manual);
g_hash_table_destroy(sysenv.image_table);
module_free();

#ifndef _WIN32
host_free_all();
#endif
}

/**********************************/
/* read the init file if possible */
/**********************************/
gint read_gdisrc(void)
{
gint i, element_flag, num_tokens;
gchar line[LINELEN], **buff;
gdouble version;
FILE *fp;

/* attempt to open */
fp = fopen(sysenv.init_file, "r");

/* check for an old/bad gdisrc */
if (fp)
  {
/* scan the 1st line */
  buff = get_tokenized_line(fp, &num_tokens);
  if (g_ascii_strncasecmp(*buff,"gdis",4) == 0)
    {
/* test for right version */
    if (num_tokens < 2)
      return(1);
    version = str_to_float(*(buff+1));
    sscanf(*(buff+1),"%lf",&version);
/* 0.75 was when the new element scheme was implemented */
    if (version < 0.75)
      return(1);
    }
  else
    return(1);
  g_strfreev(buff);
  }
else
  return(1);

/* read until EOF, or (failsafe) reached elements[] array allocation */
element_flag=0;
i=0;
for (;;)
  {
/* NB: we need line[] to get the OpenGL font */
  if (fgetline(fp, line))
    break;
  buff = tokenize(line, &num_tokens);
  if (!buff)
    break;

/* decide what to read */
  if (g_ascii_strncasecmp("size",*buff,4) == 0)
    {
    if (num_tokens > 2)
      {
      sysenv.width = str_to_float(*(buff+1));
      sysenv.height = str_to_float(*(buff+2));
      if (sysenv.width < MIN_WIDTH)
        sysenv.width = MIN_WIDTH;
      if (sysenv.height < MIN_HEIGHT)
        sysenv.height = MIN_HEIGHT;
      }
    } 

/* model pane width */
  if (g_ascii_strncasecmp("pane",*buff,4) == 0)
    {
    if (num_tokens > 2)
      {
      sysenv.tree_width = str_to_float(*(buff+1));
      sysenv.tray_height = str_to_float(*(buff+2));
      }
    }

/* model tree/property divider */
  if (g_ascii_strncasecmp("divide",*buff,6) == 0)
    {
    if (num_tokens > 1)
      sysenv.tree_divider = (gint) str_to_float(*(buff+1));
    }

/* povray executable */
  if (g_ascii_strncasecmp("povray_p", *buff, 8) == 0)
    sysenv.povray_path = parse_strip_newline(&line[12]);

  if (g_ascii_strncasecmp("povray_e", *buff, 8) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.povray_exe);
      sysenv.povray_exe = g_strdup(*(buff+1));
      }
    }

/* animation creation tool */
  if (g_ascii_strncasecmp("convert_p", *buff, 9) == 0)
    sysenv.convert_path = parse_strip_newline(&line[13]);

  if (g_ascii_strncasecmp("convert_e", *buff, 9) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.convert_exe);
      sysenv.convert_exe = g_strdup(*(buff+1));
      }
    }

/* image viewing */
  if (g_ascii_strncasecmp("viewer_p", *buff, 8) == 0)
    sysenv.viewer_path = parse_strip_newline(&line[12]);

  if (g_ascii_strncasecmp("viewer_e", *buff, 8) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.viewer_exe);
      sysenv.viewer_exe = g_strdup(*(buff+1));
      }
    }

/* GULP */
  if (g_ascii_strncasecmp("gulp_p", *buff, 6) == 0)
    sysenv.gulp_path = parse_strip_newline(&line[10]);

  if (g_ascii_strncasecmp("gulp_e", *buff, 6) == 0)
    {
    if (num_tokens > 1)
      {
      g_free(sysenv.gulp_exe);
      sysenv.gulp_exe = g_strdup(*(buff+1));
      }
    }
/* OpenGL drawing font */
  if (g_ascii_strncasecmp("gl_font",*buff,7) == 0)
    if (num_tokens > 1)
      strcpy(sysenv.gl_fontname, g_strstrip(&line[8]));

/* model tree box */
  if (g_ascii_strncasecmp("mtb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.mtb_on = (gint) str_to_float(*(buff+1));

/* model properties box */
  if (g_ascii_strncasecmp("mpb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.mpb_on = (gint) str_to_float(*(buff+1));

/* model symmetry box */
  if (g_ascii_strncasecmp("msb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.msb_on = (gint) str_to_float(*(buff+1));

/* atom properties box */
  if (g_ascii_strncasecmp("apb",*buff,3) == 0)
    if (num_tokens > 1)
      sysenv.apb_on = (gint) str_to_float(*(buff+1));

/* halo type */
  if (g_ascii_strncasecmp("halo",*buff,4) == 0)
    if (num_tokens > 1)
      sysenv.render.halos = (gint) str_to_float(*(buff+1));

/* low quality rotation */
  if (g_ascii_strncasecmp("fast",*buff,4) == 0)
    if (num_tokens > 1)
      sysenv.render.fast_rotation = (gint) str_to_float(*(buff+1));

  if (g_ascii_strncasecmp(*buff, "colour_bg", 9) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.bg_colour[0] = str_to_float(*(buff+1));
      sysenv.render.bg_colour[1] = str_to_float(*(buff+2));
      sysenv.render.bg_colour[2] = str_to_float(*(buff+3));
      }
    }
  if (g_ascii_strncasecmp(*buff, "colour_morph", 11) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.morph_colour[0] = str_to_float(*(buff+1));
      sysenv.render.morph_colour[1] = str_to_float(*(buff+2));
      sysenv.render.morph_colour[2] = str_to_float(*(buff+3));
      }
    }
  if (g_ascii_strncasecmp("colour_rsurf", *buff, 12) == 0)
    {
    if (num_tokens > 3)
      {
      sysenv.render.rsurf_colour[0] = str_to_float(*(buff+1));
      sysenv.render.rsurf_colour[1] = str_to_float(*(buff+2));
      sysenv.render.rsurf_colour[2] = str_to_float(*(buff+3));
      }
    }

/* cleanup */
  g_strfreev(buff);
  }

/* parse for elements */
rewind(fp);
read_elem_data(fp, MODIFIED);

return(0);
}

/*********************************************/
/* write setup & elements to the gdisrc file */
/*********************************************/
gint write_gdisrc(void)
{
FILE *fp;

fp = fopen(sysenv.init_file,"w");
if (!fp)
  {
  printf("Error: unable to create %s\n", sysenv.init_file);
  return(1);
  }

/* save final dimensions */
if (sysenv.mpane)
  if (GTK_IS_WIDGET(sysenv.mpane))
    sysenv.tree_width = GTK_WIDGET(sysenv.mpane)->allocation.width;
if (sysenv.tpane)
  if (GTK_IS_WIDGET(sysenv.tpane))
    sysenv.tray_height = GTK_WIDGET(sysenv.tpane)->allocation.height;

fprintf(fp,"gdis %f\n", VERSION);
fprintf(fp,"canvas %d\n", sysenv.canvas);
fprintf(fp,"size %d %d\n", sysenv.width,sysenv.height);
fprintf(fp,"pane %d %d\n", sysenv.tree_width, sysenv.tray_height);
fprintf(fp,"divider %d\n", sysenv.tree_divider);
fprintf(fp,"gl_font %s\n", sysenv.gl_fontname);
fprintf(fp,"mtb %d\n", sysenv.mtb_on);
fprintf(fp,"mpb %d\n", sysenv.mpb_on);
fprintf(fp,"msb %d\n", sysenv.msb_on);
fprintf(fp,"apb %d\n", sysenv.apb_on);
fprintf(fp,"halos %d\n", sysenv.render.halos);
fprintf(fp,"fast %d\n", sysenv.render.fast_rotation);

fprintf(fp,"colour_bg %f %f %f\n", sysenv.render.bg_colour[0],
                                   sysenv.render.bg_colour[1],
                                   sysenv.render.bg_colour[2]);
fprintf(fp,"colour_morph %f %f %f\n", sysenv.render.morph_colour[0],
                                      sysenv.render.morph_colour[1],
                                      sysenv.render.morph_colour[2]);
fprintf(fp,"colour_rsurf %f %f %f\n", sysenv.render.rsurf_colour[0],
                                      sysenv.render.rsurf_colour[1],
                                      sysenv.render.rsurf_colour[2]);

if (sysenv.babel_path)
  fprintf(fp,"babel_path %s\n", sysenv.babel_path);
if (sysenv.convert_path)
  fprintf(fp,"convert_path %s\n", sysenv.convert_path);
if (sysenv.gulp_path)
  fprintf(fp,"gulp_path %s\n", sysenv.gulp_path);
if (sysenv.gamess_path)
  fprintf(fp,"gamess_path %s\n", sysenv.gamess_path);
if (sysenv.povray_path)
  fprintf(fp,"povray_path %s\n", sysenv.povray_path);
if (sysenv.viewer_path)
  fprintf(fp,"viewer_path %s\n", sysenv.viewer_path);


fprintf(fp,"babel_exe %s\n", sysenv.babel_exe);
fprintf(fp,"convert_exe %s\n", sysenv.convert_exe);
fprintf(fp,"gulp_exe %s\n", sysenv.gulp_exe);
fprintf(fp,"gamess_exe %s\n", sysenv.gamess_exe);
fprintf(fp,"povray_exe %s\n", sysenv.povray_exe);
fprintf(fp,"viewer_exe %s\n", sysenv.viewer_exe);



/* write the non-default element data */
write_elem_data(fp);

fclose(fp);
return(0);
}

/**************/
/* main setup */
/**************/
#define DEBUG_SYS_INIT 0
void sys_init(gint argc, gchar *argv[])
{
gchar *temp;
const gchar *ctemp;
struct light_pak *light;
FILE *fp;

/* top level structure initialization */
sysenv.max_threads = 1;
sysenv.task_list = NULL;
sysenv.host_list = NULL;
sysenv.dialog_list = NULL;
sysenv.glconfig = NULL;
sysenv.stereo = FALSE;
sysenv.stereo_windowed = FALSE;
sysenv.stereo_fullscreen = FALSE;
sysenv.canvas_list = NULL;
sysenv.canvas_rows = 1;
sysenv.canvas_cols = 1;
sysenv.width = START_WIDTH;
sysenv.height = START_HEIGHT;
sysenv.snapshot = FALSE;
sysenv.snapshot_filename = NULL;
sysenv.tree_width = START_WIDTH/4;
sysenv.tree_divider = -1;
sysenv.tray_height = 60;
sysenv.mpane = NULL;
sysenv.tpane = NULL;
sysenv.fps = 40;
sysenv.moving = FALSE;
sysenv.roving = FALSE;
sysenv.refresh_canvas = FALSE;
sysenv.refresh_dialog = FALSE;
sysenv.refresh_tree = FALSE;
sysenv.refresh_properties = FALSE;
sysenv.refresh_text = FALSE;
sysenv.mtb_on = TRUE;
sysenv.mpb_on = TRUE;
sysenv.msb_on = TRUE;
sysenv.apb_on = TRUE;
sysenv.lmb_on = FALSE;
sysenv.pib_on = FALSE;
/* default to single model display */
sysenv.mal = NULL;
sysenv.displayed[0] = -1;
sysenv.active_model = NULL;
sysenv.canvas = TRUE;
sysenv.select_mode = CORE;
sysenv.calc_pbonds = TRUE;

/* default masks */
sysenv.file_type = DATA;
sysenv.babel_type = AUTO;
sysenv.num_elements = sizeof(elements) / sizeof(struct elem_pak);
sysenv.elements = NULL;

/* rendering setup */
sysenv.render.width = 600;
sysenv.render.height = 600;
sysenv.render.vp_dist = 5.0;
sysenv.render.zone_size = 10.0;
sysenv.render.show_energy = FALSE;

/* stereo defaults */
sysenv.render.stereo_quadbuffer = FALSE;
sysenv.render.stereo_use_frustum = TRUE;
sysenv.render.stereo_eye_offset = 1.0;
sysenv.render.stereo_parallax = 1.0;
sysenv.render.stereo_left = TRUE;
sysenv.render.stereo_right = TRUE;

/* default light */
light = g_malloc(sizeof(struct light_pak));
light->type = DIRECTIONAL;
VEC3SET(light->x, -1.0, 1.0, -1.0);
VEC3SET(light->colour, 1.0, 1.0, 1.0);
light->ambient = 0.2;
light->diffuse = 0.8;
light->specular = 0.7;
/*** JJM new values for better (IMHO) lighting ***/
light->diffuse = 0.5;
//light->specular = 0.3;
sysenv.render.light_list = g_slist_append(sysenv.render.light_list, light);
sysenv.render.num_lights = 1;

sysenv.render.perspective = FALSE;
sysenv.render.antialias = FALSE;
sysenv.render.wire_show_hidden = FALSE;
sysenv.render.fog = FALSE;
sysenv.render.wire_model = FALSE;
sysenv.render.wire_surface = FALSE;
sysenv.render.shadowless = FALSE;
sysenv.render.animate = FALSE;
sysenv.render.no_povray_exec = FALSE;
sysenv.render.no_keep_tempfiles = TRUE;
sysenv.render.animate_type = ANIM_GIF;
sysenv.render.animate_file = g_strdup("movie");
sysenv.render.atype = FALSE;
sysenv.render.axes = FALSE;
sysenv.render.morph_finish = g_strdup("F_Glass4");
sysenv.render.ref_index = 2.5;
sysenv.render.delay = 20.0;
sysenv.render.mpeg_quality = 95.0;
sysenv.render.sphere_quality = 3.0;
sysenv.render.cylinder_quality = 9.0;
sysenv.render.ribbon_quality = 10.0;
sysenv.render.ms_grid_size = 0.5;
sysenv.render.auto_quality = FALSE;
sysenv.render.fast_rotation = FALSE;
sysenv.render.halos = FALSE;
sysenv.render.ambience = 0.2;
sysenv.render.diffuse = 0.9;
sysenv.render.specular = 0.9;
sysenv.render.transmit = 1.0;
sysenv.render.ghost_opacity = 0.5;
sysenv.render.ball_rad = 0.3;
sysenv.render.stick_rad = 0.1;
sysenv.render.stick_thickness = GTKGL_LINE_WIDTH;
sysenv.render.line_thickness = GTKGL_LINE_WIDTH;
sysenv.render.frame_thickness = GTKGL_LINE_WIDTH;
sysenv.render.geom_line_width = 3.0;
sysenv.render.cpk_scale = 1.0;
sysenv.render.fog_density = 0.35;
sysenv.render.fog_start = 0.5;
sysenv.render.ribbon_curvature = 0.5;
sysenv.render.ribbon_thickness = 1.0;
sysenv.render.phonon_scaling = 4.0;
sysenv.render.ahl_strength = 0.7;
sysenv.render.ahl_size = 20;
sysenv.render.shl_strength = 0.8;
sysenv.render.shl_size = 20;
VEC3SET(sysenv.render.fg_colour, 1.0, 1.0, 1.0);
VEC3SET(sysenv.render.bg_colour, 0.0, 0.0, 0.0);
VEC3SET(sysenv.render.morph_colour, 0.1, 0.1, 0.8);
VEC3SET(sysenv.render.rsurf_colour, 0.0, 0.3, 0.6);
VEC3SET(sysenv.render.label_colour, 1.0, 1.0, 0.0);
VEC3SET(sysenv.render.title_colour, 0.0, 1.0, 1.0);
VEC3SET(sysenv.render.ribbon_colour, 0.0, 0.0, 1.0);

/* setup directory and file pointers */
sysenv.cwd = g_get_current_dir();
const gchar *envdir = g_getenv("GDIS_START_DIR");
if (envdir)
  sysenv.cwd = (gchar *) envdir;

#if DEBUG_SYS_INIT
printf("cwd: [%s]\n", sysenv.cwd);
#endif

/* generate element file full pathname */
/* sometimes this returns the program name, and sometimes it doesn't */
temp = g_find_program_in_path(argv[0]);
/* remove program name (if attached) */
if (g_file_test(temp, G_FILE_TEST_IS_DIR))
  sysenv.gdis_path = temp;
else
  {
  sysenv.gdis_path = g_path_get_dirname(temp);
  g_free(temp);
  }

#if DEBUG_SYS_INIT
printf("%s path: [%s]\n", argv[0], sysenv.gdis_path);
#endif

if (sysenv.gdis_path)
  sysenv.elem_file = g_build_filename(sysenv.gdis_path, ELEM_FILE, NULL);
else
  {
  printf("WARNING: gdis directory not found.\n");
  sysenv.elem_file = g_build_filename(sysenv.cwd, ELEM_FILE, NULL);
  }

/* generate gdisrc full pathname */
ctemp = g_get_home_dir();
if (ctemp)
  sysenv.init_file = g_build_filename(ctemp, INIT_FILE, NULL);
else
  {
  printf("WARNING: home directory not found.\n");
  if (sysenv.gdis_path)
    sysenv.init_file = g_build_filename(sysenv.gdis_path, INIT_FILE, NULL);
  else
    {
    if (sysenv.cwd)
      sysenv.init_file = g_build_filename(sysenv.cwd, INIT_FILE, NULL);
    else
      {
      printf("FATAL: current directory not found.\n");
      exit(-1);
      }
    }
  }

/* defaults */
#if _WIN32
sysenv.babel_exe = g_strdup("babel.exe");
sysenv.povray_exe = g_strdup("povray.exe");
sysenv.convert_exe = g_strdup("convert.exe");
sysenv.viewer_exe = g_strdup("display.exe");
sysenv.gulp_exe = g_strdup("gulp.exe");
sysenv.gamess_exe = g_strdup("wingamess");
#else
sysenv.babel_exe = g_strdup("babel");
sysenv.povray_exe = g_strdup("povray");
sysenv.convert_exe = g_strdup("convert");
sysenv.viewer_exe = g_strdup("display");
sysenv.gulp_exe = g_strdup("gulp");
sysenv.gamess_exe = g_strdup("run_gms_for_gdis");
#endif


strcpy(sysenv.gl_fontname, GL_FONT);

/* atomic scattering factor coefficients */
sysenv.sfc_table = g_hash_table_new_full(&g_str_hash,
                                         &hash_strcmp,
                                         &g_free,
                                         &free_sfc);

/* IMPORTANT this must be done before gdisrc is parsed as */
/* setup the element data before any possible modifications */
/* can be read in from the gdisrc file */
printf("scanning: %-50s ", sysenv.elem_file);
fp = fopen(sysenv.elem_file, "rt");
if (fp)
  {
  read_elem_data(fp, DEFAULT);
  fclose(fp);
  printf("[ok]\n");
  }
else
  {
/* missing default element info is fatal */
  printf("[not found]\n");
  exit(1);
  }

/* if failed to read gdisrc - write a new one */
/* TODO - if overwrite old version, save/rename first? */
sysenv.babel_path = NULL;
sysenv.convert_path = NULL;
sysenv.gulp_path = NULL;
sysenv.gamess_path = NULL;
sysenv.povray_path = NULL;
sysenv.viewer_path = NULL;
if (read_gdisrc())
  {
  printf("creating: %-50s ", sysenv.init_file);
  if (write_gdisrc())
    printf("[error]\n");
  else
    printf("[ok]\n");
  }
else
  printf("scanning: %-50s [ok]\n", sysenv.init_file);

/* get executable paths */
if (!sysenv.babel_path)
  sysenv.babel_path = g_find_program_in_path(sysenv.babel_exe);
if (!sysenv.convert_path)
  sysenv.convert_path = g_find_program_in_path(sysenv.convert_exe);
if (!sysenv.gulp_path)
  sysenv.gulp_path = g_find_program_in_path(sysenv.gulp_exe);
if (!sysenv.povray_path)
  sysenv.povray_path = g_find_program_in_path(sysenv.povray_exe);

/* display program */
if (!sysenv.viewer_path)
  sysenv.viewer_path = g_find_program_in_path(sysenv.viewer_exe);

/* alternative (FIXME - starting to get ugly - redo) */
if (!sysenv.viewer_path)
  sysenv.viewer_path = g_find_program_in_path("eog");

/* HACK FIX - GAMESS path contains the path, not the path + file */
temp = g_find_program_in_path(sysenv.gamess_exe);

if (temp)
  sysenv.gamess_path = g_path_get_dirname(temp);
else
#if __WIN32
  sysenv.gamess_path = g_strdup("c:\\wingamess");
#else
  sysenv.gamess_path = NULL;
#endif

g_free(temp);

/* tree init */
sysenv.num_trees = 0;
/* copied selection */
sysenv.select_source = NULL;
sysenv.module_list = NULL;
sysenv.projects = NULL;
sysenv.manual = NULL;
sysenv.image_table = NULL;
sysenv.surfaces = NULL;

/* module */
file_init();
#ifdef WITH_GUI
help_init();
#endif
init_math();
library_init();
}

/********/
/* MAIN */
/********/
gint main(int argc, char *argv[])
{
/* NBL read this 1st as it affects canvas type, but command */
/* line arguments should be able to overide */
sys_init(argc, argv);
module_setup();
sysenv.write_gdisrc = FALSE;

/* command line? */
if (argc > 1)
  if (g_ascii_strncasecmp(argv[1],"-c",3) == 0)
    sysenv.canvas = FALSE;


#ifdef WITH_GUI
/* set up main window and event handlers (or start non-gui mode) */
if (sysenv.canvas)
  {
  gint i;

/* initialize threads and set up the queue */
  task_queue_init();

/* CURRENT - absorb these into gui_init() */
/*
  gtk_init(&argc, &argv);
  gtk_gl_init(&argc, &argv);
*/

/* main interface */
  gui_init(argc, argv);

/* task update timer */
/* TODO - put this in gui_widget_handler? */
  g_timeout_add(1000, (GSourceFunc) &update_task_info, NULL);

/* screen redraw timer */
  g_timeout_add(25, (GSourceFunc) &gui_canvas_handler, NULL);

/* CURRENT - gui updates done through a timer */
/* this was done so that threads can safely schedule a gui update */
  g_timeout_add(200, (GSourceFunc) &gui_widget_handler, NULL);

/* process arguments as files to load */
  for (i=1 ; i<argc ; i++)
    file_load(argv[i], NULL);

/* thread-safe handling */
  gdk_threads_enter();
  gtk_main();
  gdk_threads_leave();
  }
else
#else
sysenv.canvas = FALSE;
#endif
  command_main_loop(argc, argv);

return(0);
}

/* CURRENT */
/* routines that are not cleanly seperable when we build with no GUI */
#ifndef WITH_GUI
void gui_text_show(gint type, gchar *msg)
{
printf("%s", msg);
}

void gui_refresh(gint dummy)
{
}

void tree_select_delete(void)
{
}
void tree_select_active(void)
{
}
void tree_select_model(struct model_pak *m)
{
}
void tree_model_add(struct model_pak *m)
{
}
void tree_model_refresh(struct model_pak *m)
{
}

void dialog_destroy_type(gint dummy)
{
}

gpointer graph_new(const gchar *dummy, struct model_pak *m)
{
return(NULL);
}

void graph_add_data(gint a, gdouble *b, gdouble c, gdouble d, gpointer e)
{
}

void graph_set_yticks(gint a, gint b, gpointer c)
{
}

void graph_free_list(struct model_pak *m)
{
}

void meas_graft_model(struct model_pak *m)
{
}

void meas_prune_model(struct model_pak *m)
{
}

/* FIXME - this should actually be replaced by gui_refresh() */
void redraw_canvas(gint dummy)
{
}
#endif
