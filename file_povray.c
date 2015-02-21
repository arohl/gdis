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
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

#include "gdis.h"
#include "coords.h"
#include "file.h"
#include "matrix.h"
#include "quaternion.h"
#include "numeric.h"
#include "morph.h"
#include "opengl.h"
#include "render.h"
#include "select.h"
#include "spatial.h"
#include "interface.h"
#include "colourlib.h"

extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];

/* global render parameters */
gdouble p2a;

/***************************/
/* Setup POVray parameters */
/***************************/
void povray_hdr(FILE *fp, struct model_pak *data)
{
gdouble xvec, yvec, amb, pos[3], colour[3];
gdouble x[3], o[3], v[3], e[3];
GSList *list;
struct light_pak *light;
struct camera_pak *camera;

g_assert(data != NULL);
g_assert(data->camera != NULL);

fprintf(fp,"#include \"colors.inc\" \n");
fprintf(fp,"#include \"finish.inc\" \n");
fprintf(fp,"#include \"glass.inc\" \n");
fprintf(fp,"#include \"metals.inc\" \n");
fprintf(fp,"#include \"textures.inc\" \n");

/* background colour (except for glass morphologies) */
fprintf(fp,"background { color rgb<%f,%f,%f0> }\n", sysenv.render.bg_colour[0],
                        sysenv.render.bg_colour[1], sysenv.render.bg_colour[2]);  

/* pixel to angstrom conversion, with yet another magic number... */
p2a = 0.565 * (gdouble) sysenv.render.width / data->rmax;

/* preserve model aspect ratio for the given image size */
xvec = yvec = 2.0*sysenv.rsize;
if (sysenv.render.width > sysenv.render.height)
  xvec *= sysenv.render.width/sysenv.render.height;
if (sysenv.render.height > sysenv.render.width)
  yvec *= sysenv.render.height/sysenv.render.width;

/* compute camera position and orientation */
camera = data->camera;
ARR3SET(x, camera->x);
ARR3SET(o, camera->o);
ARR3SET(v, camera->v);

switch (camera->mode)
  {
  case FREE:
    break;

  default:
  case LOCKED:
    quat_rotate(x, camera->q);
    quat_rotate(o, camera->q);
    quat_rotate(v, camera->q);
    break;
  }
/* convert viewing vector to a location */
ARR3ADD(v, x);

/* camera zoom */
xvec *= camera->zoom;
yvec *= camera->zoom;

/* NEW - enable movies of left/right eye to be produced */
if (sysenv.stereo)
  {
/* get axis for eye translation (view x up vector) */
  crossprod(e, v, o);
  normalize(e, 3);

/* the old 2% rule ... */
  VEC3MUL(e, 0.02 * sysenv.rsize); 

/* default is left eye only */
  if (sysenv.render.stereo_right)
    {
    ARR3ADD(x, e);
    ARR3ADD(v, e);
    }
  else
    {
    ARR3SUB(x, e);
    ARR3SUB(v, e);
    }
  }

/* sky is the orientation vector */
/* right and up give the perspective */
if (camera->perspective)
  {
  fprintf(fp,"camera { location <%f,%f,%f>\n", x[0], x[1], x[2]);
  fprintf(fp,"    sky <%f,%f,%f>\n", o[0], o[1], o[2]);
  fprintf(fp,"    right <%f,0,0> up <%f,0,0>\n", xvec, yvec);
  fprintf(fp,"    look_at <%f,%f,%f>\n", v[0], v[1], v[2]);
  fprintf(fp,"    angle %f }\n", camera->fov);
  }
else
  {
  fprintf(fp,"camera { orthographic location <%f,%f,%f>\n", x[0], x[1], x[2]);
  fprintf(fp,"    sky <%f,%f,%f>\n", o[0], o[1], o[2]);
  fprintf(fp,"    right <%f,0,0> up <%f,0,0>\n", xvec, yvec);
  fprintf(fp,"    look_at <%f,%f,%f> }\n", v[0], v[1], v[2]);
  }

/* create light sources */
for (list=sysenv.render.light_list ; list ; list=g_slist_next(list))
  {
  light = list->data;
  ARR3SET(pos, light->x);

/* OpenGL -> POVRay axes */
  pos[0] *= -1.0;
  pos[1] *= -1.0;
  pos[2] *= -1.0;

  quat_rotate(pos, camera->q);

  switch (light->type)
    {
    case POSITIONAL:
      fprintf(fp,"light_source\n  {\n <%f,%f,%f>\n", pos[0], pos[1], pos[2]);
      break;

    case DIRECTIONAL:
/* move away far enough so the rays are ~ // */
      VEC3MUL(pos, 100.0*data->rmax);
      fprintf(fp,"light_source\n  {\n <%f,%f,%f>\n", pos[0], pos[1], pos[2]);
      break;

    default:
      continue;
    }

/* simulate OpenGL style lights */
  ARR3SET(colour, light->colour);
  VEC3MUL(colour, light->specular);

  if (sysenv.render.shadowless)
    fprintf(fp,"  color rgb<%f,%f,%f> shadowless }\n", colour[0], colour[1], colour[2]);
  else
    {
/* old style lighting */
/*
    fprintf(fp,"  color rgb<%f,%f,%f> }\n", colour[0], colour[1], colour[2]);
*/
    fprintf (fp,"color White\n");
    fprintf (fp,"   area_light <5, 0, 0,>, <0, 0, 5>, 5, 5\n");
    fprintf (fp,"   adaptive 1\n   jitter\n}\n");
    }
  }

/* fill-light to bring out the shadows */
fprintf(fp,"light_source{<%f,%f,%f> color Gray80 shadowless}\n", pos[0], pos[1], pos[2]);

/* morph is too dark with just the above, sky_sphere is *nice* */
/* TODO - choice of colour eg white/grey/light blue */
/* NB: white is a bit too bright (even with 0 ambience) */
if (data->id == MORPH)
  {
  if (!sysenv.render.wire_surface && !sysenv.render.shadowless)
    {
    fprintf(fp,"sky_sphere { pigment {gradient y  color_map "
               "{[0, 1 color Gray20 color White]} rotate x*45}}\n");
    }
  }

/* POVRay is a bit darker than OpenGL */
amb = 20.0*sysenv.render.ambience;

fprintf(fp,"global_settings { ambient_light rgb<%f, %f, %f> assumed_gamma 2.2}\n",amb,amb,amb);
}

/* TODO -rename these something sensible */
/*********************/
/* forward algorithm */
/*********************/
void povray_rgb_decode(gint n, gdouble *c)
{
c[0] = ((n >> 10) & 31)/31.0;
c[1] = ((n >> 5) & 31)/31.0;
c[2] = (n & 31)/31.0;
}

/**********************/
/* backward algorithm */
/**********************/
gint povray_rgb_encode(gdouble *c)
{
gint n;

n = nearest_int(c[2]*31);
n += nearest_int(c[1]*31) << 5;
n += nearest_int(c[0]*31) << 10;

return(n);
}

/********************************************/
/* write a set of pre-defined textures to   */
/* a POV-Ray file. These textures describe  */
/* 256 HSV colours for linearly mapping     */
/* a property from red to blue              */
/********************************************/
void write_povray_colour_textures(FILE *povfile, struct model_pak *data, int colour_style)
{
int i;
float hue, red,green,blue;
  
switch (colour_style)
  {
  case HSV:
    for (i=0; i<360; i++)
      {
      hue = (float)i;
      hsv2rgb(hue,1.0,1.0,&red,&green,&blue);
      fprintf(povfile,"#declare HSV%d = texture{pigment{color rgbt<%.2f,%.2f,%.2f,%.1f>}",i,red,green,blue,1.0-sysenv.render.transmit);
      fprintf(povfile," finish { Phong_Shiny } }\n");
      }
    break;

/* FIXME - we only need to do this texture declaration for meshes (ie smooth tri spatials) */
  case REDWHITEBLUE:
    for (i=0 ; i<32768 ; i++)
      {
      gdouble c[3];

      fprintf(povfile,"#declare RGB_%d = texture{pigment{color rgb ", i);
      povray_rgb_decode(i, c);
      fprintf(povfile," <%f,%f,%f> } finish { Phong_Shiny } }\n", c[0], c[1], c[2]);

/* debug */
/*
{
gint j;

printf("[%d] : <%f,%f,%f> :", i, c[0], c[1], c[2]);
j = povray_rgb_encode(c);
printf("[%d]\n", j);
}
*/

      }
    break;
  }
}

/****************************/
/* make a POVRAY input file */
/****************************/
/* TODO - eliminate redundancy between this and the opengl code */
/* TODO - put this stuff in file_povray.c */

/* Note by JJM. Need to minimise unnecessary blankspaces in the POV-Ray file - it can get very big! */

#define DEBUG_MKPOV 0
gint write_povray(gchar *povfile, struct model_pak *data)
{
gint n, m, i, j;
gint r, g, b;
gint style;				/* for selecting the colour style of iso-surfaces */
gdouble rad, scale, len;
gdouble x1, y1, z1, x2, y2, z2;
gdouble xmin, ymin, zmin, xmax, ymax, zmax;
gdouble rf, gf, bf;
gdouble vec[3], vec1[3], vec2[3];
gdouble v1[4], v2[4], v3[3], n1[3], n2[3], n3[3];
gdouble ctrl[8][3], mat4[16];
gint colour_index[3];
GSList *list, *ilist, *plist, *rlist, *pipe_list[4];
struct core_pak *core1;
struct object_pak *object;
struct ribbon_pak *ribbon;
struct spatial_pak *spatial;
struct pipe_pak *pipe;
struct image_pak *image;
struct vec_pak *p1, *p2, *p3;
struct elem_pak elem;
FILE *fp;

/* checks */
g_assert(data != NULL);

/* open file & init */
fp = fopen(povfile, "w");
if (!fp)
  {
  printf("Error, render(): can't open %s!\n", povfile);
  return(1);
  }

/* strcpy(povfile,models[model].filename); */
/* file_extension -> .pov */

/* setup */
/* JJM temporary solution for picking the colour-style of surfaces */
//style = HSV;
style = REDWHITEBLUE;

povray_hdr(fp,data);

write_povray_colour_textures(fp,data,style);

/* limits - for intelligent axes placement */
xmin = ymin = zmin = 99999999.9;
xmax = ymax = zmax = -99999999.9;
/* pixel to coord scaling factor */
scale = sysenv.subscale * 0.5 / data->rmax;

/* FIXME - can we just sub in data->display_lattice??? */
/* precalc matrix products */
ARR3SET(&mat4[0], &data->latmat[0]);
ARR3SET(&mat4[4], &data->latmat[3]);
ARR3SET(&mat4[8], &data->latmat[6]);
ARR3SET(v1, data->centroid);
vecmat(data->latmat, v1);
mat4[3] = -v1[0];
mat4[7] = -v1[1];
mat4[11] = -v1[2];
VEC4SET(&mat4[12], 0.0, 0.0, 0.0, 1.0);

/* do atoms */ 
#if DEBUG_MKPOV
printf("Doing atoms...\n");
#endif
/* enumerate periodic images */
plist=NULL;
do
  {
  if (plist)
    {
/* image */
    image = plist->data;
    ARR3SET(vec, image->rx);
    plist = g_slist_next(plist);
    }
  else
    {
/* original */
    VEC3SET(vec, 0.0, 0.0, 0.0);
    plist = data->images;
    }

/* only if same type */
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core1 = list->data;
  if (core1->status & (DELETED | HIDDEN))
    continue;

  ARR3SET(vec1, core1->rx);
  ARR3ADD(vec1, vec);

/* get current colour */
    r = core1->colour[0];
    g = core1->colour[1];
    b = core1->colour[2];
/* convert to povray rgb format */
    rf = (gdouble) (r) / 65535.0;
    gf = (gdouble) (g) / 65535.0;
    bf = (gdouble) (b) / 65535.0;
    switch(core1->render_mode)
      {
/* rounded end are now draw in the bond section */
      case LIQUORICE:
      case STICK:
        if (core1->bonds)
          break;
      case BALL_STICK:
        fprintf(fp,"sphere { <%f, %f, %f>, %f ",vec1[0],vec1[1],vec1[2],
                                             sysenv.render.ball_rad);
/* TODO - can we adjust this to get a wire frame sphere? */
        fprintf(fp,"texture{pigment{color rgb<%f,%f,%f>",rf,gf,bf);  

        if (core1->ghost)
          fprintf(fp, " transmit 0.6}\n");
        else
          fprintf(fp, "}\n");

		fprintf(fp, " finish{Phong_Shiny}}}\n");
/*        fprintf(fp,"  finish{phong %f phong_size %d }}}\n",
                    sysenv.render.ahl_strength, (gint) sysenv.render.ahl_size); */
/*
        fprintf(fp,"  finish{specular 0.6 diffuse 0.8 }}}\n");
*/

        break;

      case CPK:
/* FIXME - calling get_elem_data() all the time is inefficient */
        get_elem_data(core1->atom_code, &elem, data);
        rad = elem.vdw;
        rad *= sysenv.render.cpk_scale;
        fprintf(fp,"sphere { <%f, %f, %f>, %f ",vec1[0],vec1[1],vec1[2],rad);
        fprintf(fp,"texture{pigment{color rgb<%f,%f,%f>",rf,gf,bf);  

        if (core1->ghost)
          fprintf(fp, " transmit 0.6}\n");
        else
          fprintf(fp, "}\n");

		fprintf(fp, " finish{Phong_Shiny}}}\n");
/*        fprintf(fp,"  finish{phong %f phong_size %d }}}\n",
                    sysenv.render.ahl_strength, (gint) sysenv.render.ahl_size); */
/*
        fprintf(fp,"  finish{specular 0.6 diffuse 0.8 }}}\n");
*/
        break;
      }
  }
  }
while (plist);

#if DEBUG_MKPOV
printf("Doing bonds...\n");
#endif

/* calculate all bonds */
render_make_pipes(pipe_list, data);

/* enumerate the supplied pipes (half bonds) */
for (i=0 ; i<4 ; i++)
  {
  for (list=pipe_list[i] ; list ; list=g_slist_next(list))
    {
    pipe = list->data;

/* original + image iteration */
    ilist = NULL;
    do
      {
/* original */
      ARR3SET(v1, pipe->v1);
      ARR3SET(v2, pipe->v2);
      if (ilist)
        {
        image = ilist->data;
/* image */
        ARR3ADD(v1, image->rx);
        ARR3ADD(v2, image->rx);
        ilist = g_slist_next(ilist);
        }
      else
        ilist = data->images;

/* bond type */
      switch(i)
        {
/* normal */
        case 0:
          fprintf(fp,"cylinder { <%f,%f,%f>,\n<%f,%f,%f>, %f\n",
                      v1[0],v1[1],v1[2],v2[0],v2[1],v2[2], sysenv.render.stick_rad);
          fprintf(fp,"open texture{pigment{color rgb<%f,%f,%f>}\n",
                      pipe->colour[0], pipe->colour[1], pipe->colour[2]);
		  fprintf(fp, " finish{Phong_Shiny}}}\n");
/*          fprintf(fp,"  finish{phong %f phong_size %d }}}\n",
                      sysenv.render.ahl_strength, (gint) sysenv.render.ahl_size); */
          break;

/* ghosts */
        case 1:
          fprintf(fp,"cylinder { <%f,%f,%f>,\n<%f,%f,%f>, %f\n",
                      v1[0],v1[1],v1[2],v2[0],v2[1],v2[2], sysenv.render.stick_rad);
          fprintf(fp,"open texture{pigment{color rgb<%f,%f,%f> transmit 0.6}",
                      pipe->colour[0], pipe->colour[1], pipe->colour[2]);
		  fprintf(fp, " finish{Phong_Shiny}}}\n");
/*          fprintf(fp,"  finish{phong %f phong_size %d }}}\n",
                      sysenv.render.ahl_strength, (gint) sysenv.render.ahl_size); */
          break;

/* line */
        case 3:
/* FIXME - need a better way to set "line" thicknesses */
          fprintf(fp,"cylinder { <%f,%f,%f>,\n<%f,%f,%f>, %f\n",
                      v1[0],v1[1],v1[2],v2[0],v2[1],v2[2], 2.0/p2a);
          fprintf(fp,"open texture{pigment{color rgb<%f,%f,%f>}\n",
                      pipe->colour[0], pipe->colour[1], pipe->colour[2]);
		  fprintf(fp, " finish{Phong_Shiny}}}\n");
		  /*
          fprintf(fp,"  finish{phong %f phong_size %d }}}\n",
                      sysenv.render.ahl_strength, (gint) sysenv.render.ahl_size); */
          break;

/* TODO - wire frame bonds??? */
        default:
          break;
        }
      }
    while (ilist);
    }
  }

/* free all pipes */
for (i=4 ; i-- ; )
  free_slist(pipe_list[i]);

#if DEBUG_MKPOV
printf("Doing special objects...\n");
#endif

/********/
/* AXES */
/********/
/* FIXME - stuffs up under new camera stuff */
if (0)
if (data->show_axes)
  {
/* origin */
  VEC3SET(vec1, 1.0, 1.0, 0.0);
  VEC3MUL(vec1, data->rmax);
/* cope with scaling */
  vec1[2] -= 4.0*sysenv.rsize;
  for (j=0 ; j<3 ; j++)
    {
/* get end point */
    ARR3SET(vec2, data->axes[j].rx);
    ARR3ADD(vec2, vec1);
/* cope with scaling */
    vec2[2] -= 4.0*sysenv.rsize;

/* HACK for degenerate cylinders (povray spews on these) */
ARR3SET(v3, vec1);
ARR3SUB(v3, vec2);
if (VEC3MAGSQ(v3) < 0.1)
  continue;

/* draw */
    fprintf(fp,"cylinder { <%f,%f,%f>,<%f,%f,%f>, %f\n",
                vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2], 2.0/p2a);
    fprintf(fp," texture { pigment {White} } }\n");
    }
  }

/********/
/* CELL */
/********/
if (data->show_cell && data->periodic)
  {
  rad = sysenv.render.frame_thickness / p2a;

/* ends */
  for (j=0 ; j<8 ; j++)
    {
    m = 2*(j/2) + 1;
    n = 4*(j/4);

    ARR3SET(vec, data->cell[m].rx);
    x1 = vec[0];
    y1 = vec[1];
    z1 = vec[2];

    ARR3SET(vec, data->cell[n].rx);
    x2 = vec[0];
    y2 = vec[1];
    z2 = vec[2];

/* HACK for degenerate cylinders (povray spews on these) */
ARR3SET(v3, data->cell[m].rx);
ARR3SUB(v3, data->cell[n].rx);
if (VEC3MAGSQ(v3) < 0.1)
  continue;

    fprintf(fp,"cylinder { <%f,%f,%f>,<%f,%f,%f>, %f\n",
                                 x1,y1,z1,x2,y2,z2,rad);
    fprintf(fp," texture { pigment {White} } }\n");

/* HACK for degenerate cylinders (povray spews on these) */
ARR3SET(v3, data->cell[m].rx);
ARR3SUB(v3, data->cell[n+2].rx);
if (VEC3MAGSQ(v3) < 0.1)
  continue;

    ARR3SET(vec, data->cell[n+2].rx);
    x2 = vec[0];
    y2 = vec[1];
    z2 = vec[2];
    fprintf(fp,"cylinder { <%f,%f,%f>,<%f,%f,%f>, %f\n",
                                 x1,y1,z1,x2,y2,z2,rad);
    fprintf(fp," texture { pigment {White} } }\n");
    }

/* sides */
/* skip for 2D periodic models */
  if (data->periodic == 3)
    {
    m = 0;
    n = 4;
    for (j=0 ; j<4 ; j++)
      {
      ARR3SET(vec, data->cell[m].rx);
      x1 = vec[0];
      y1 = vec[1];
      z1 = vec[2];

      ARR3SET(vec, data->cell[n].rx);
      x2 = vec[0];
      y2 = vec[1];
      z2 = vec[2];

      fprintf(fp,"cylinder { <%f,%f,%f>,<%f,%f,%f>, %f\n",
                                   x1,y1,z1,x2,y2,z2,rad);
      fprintf(fp," texture { pigment {White} } }\n");
      m++;
      n++;
      }
    }
  }

/***********/
/* ribbons */
/***********/
for (list=data->ribbons ; list ; list=g_slist_next(list))
  {
  object = list->data;
  g_assert(object->type == RIBBON);

  rlist = (GSList *) object->data;
  while (rlist)
    {
    ribbon = rlist->data;

  fprintf(fp, "bicubic_patch {\n  type 1 flatness 0.001\n");
/* NB: POVRay is very slow at rendering ribbons, so halve the quality */
  fprintf(fp, "  u_steps %d  v_steps 2\n", (gint) (sysenv.render.ribbon_quality/2.0));

/* end points */
    ARR3SET(&ctrl[0][0], ribbon->r1);
    ARR3SET(&ctrl[3][0], ribbon->r2);

/* get distance between ribbon points */
    ARR3SET(vec1, ribbon->x1);
    ARR3SUB(vec1, ribbon->x2);
    len = VEC3MAG(vec1);

/* shape control points */
    ARR3SET(&ctrl[1][0], ribbon->r1);
    ARR3SET(&ctrl[2][0], ribbon->r2);

/* segment length based curvature - controls how flat it is at the cyclic group */
    ARR3SET(vec1, ribbon->o1);
    VEC3MUL(vec1, len*sysenv.render.ribbon_curvature);
    ARR3ADD(&ctrl[1][0], vec1);
    ARR3SET(vec2, ribbon->o2);
    VEC3MUL(vec2, len*sysenv.render.ribbon_curvature);
    ARR3ADD(&ctrl[2][0], vec2);

/* compute offsets for ribbon thickness */
    crossprod(vec1, ribbon->n1, ribbon->o1);
    crossprod(vec2, ribbon->n2, ribbon->o2);
    normalize(vec1, 3);
    normalize(vec2, 3);

/* thickness vectors for the two ribbon endpoints */
    VEC3MUL(vec1, 0.5*sysenv.render.ribbon_thickness);
    VEC3MUL(vec2, 0.5*sysenv.render.ribbon_thickness);

/* ensure these are pointing the same way */
    if (via(vec1, vec2, 3) > PI/2.0)
      {
      VEC3MUL(vec2, -1.0);
      }

/* init the bottom edge control points */
    ARR3SET(&ctrl[4][0], &ctrl[0][0]);
    ARR3SET(&ctrl[5][0], &ctrl[1][0]);
    ARR3SET(&ctrl[6][0], &ctrl[2][0]);
    ARR3SET(&ctrl[7][0], &ctrl[3][0]);
/* lift points to make the top edge */
    ARR3ADD(&ctrl[0][0], vec1);
    ARR3ADD(&ctrl[1][0], vec1);
    ARR3ADD(&ctrl[2][0], vec2);
    ARR3ADD(&ctrl[3][0], vec2);
/* lower points to make the bottom edge */
    ARR3SUB(&ctrl[4][0], vec1);
    ARR3SUB(&ctrl[5][0], vec1);
    ARR3SUB(&ctrl[6][0], vec2);
    ARR3SUB(&ctrl[7][0], vec2);

    fprintf(fp, "<%f, %f, %f>\n", ctrl[0][0], ctrl[0][1], ctrl[0][2]);
    for (i=1 ; i<4 ; i++)
      fprintf(fp, ", <%f, %f, %f>\n", ctrl[i][0], ctrl[i][1], ctrl[i][2]);
    for (i=0 ; i<4 ; i++)
      fprintf(fp, ", <%f, %f, %f>\n", ctrl[i][0], ctrl[i][1], ctrl[i][2]);
    for (i=4 ; i<8 ; i++)
      fprintf(fp, ", <%f, %f, %f>\n", ctrl[i][0], ctrl[i][1], ctrl[i][2]);
    for (i=4 ; i<8 ; i++)
      fprintf(fp, ", <%f, %f, %f>\n", ctrl[i][0], ctrl[i][1], ctrl[i][2]);

    fprintf(fp," texture { pigment { color rgbt<%f, %f, %f, %f> } } no_shadow }\n",
               ribbon->colour[0], ribbon->colour[1], ribbon->colour[2],
               1.0-sysenv.render.transmit);

    rlist = g_slist_next(rlist);
    }
  }

/******************/
/* spatial planes */
/******************/
/* FIXME - currently broken for morphology display */
for (list=data->spatial ; list ; list=g_slist_next(list))
  {
  spatial = list->data;

  switch(spatial->type)
    {
/* special spatials */
    case SPATIAL_VECTOR:
      break;

/* generic spatials */
    default:
/* enumerate periodic images */
      plist=NULL;
      do
        {
        if (plist)
          {
/* image */
          image = plist->data;
          ARR3SET(vec2, image->rx);
          plist = g_slist_next(plist);
          }
        else
          {
/* original */
          VEC3SET(vec2, 0.0, 0.0, 0.0);
          plist = data->images;
          }
/* enumerate vertices */
        rlist = spatial->list;
        p1 = NULL;

/* TODO - the best way of doing this would be through smooth_triangle() and mesh() */
/* note that these two seem to have to be used together ie can't just use smooth triangle */
/* by itself.  Unfortunately, a mesh is something that can have only one colour which makes */
/* it harder to implement - a rewrite would be needed. */
/* TODO - the normal calcs are redundant at present, but may be used in the future */
        switch(spatial->method)
          {
          case GL_TRIANGLES:
            fprintf(fp, "#declare Surface=mesh{\n");
            while (rlist)
              {
              fprintf(fp, "smooth_triangle{ \n");
              for (i=3 ; i-- ; )
                {
                p1 = rlist->data;
                ARR3SET(vec1, vec2);
                ARR3ADD(vec1, p1->rx);
/* can we supply normal at the vertex as well? - YES! */
                ARR3SET(vec, p1->rn);
                fprintf(fp, "<%f,%f,%f><%f,%f,%f>\n",
                             vec1[0], vec1[1], vec1[2],
                             vec[0], vec[1], vec[2]);
                colour_index[i] = povray_rgb_encode(p1->colour);
                rlist = g_slist_next(rlist);
                }
              fprintf(fp,"texture_list { ");
              fprintf(fp,"RGB_%d RGB_%d RGB_%d ",
                      colour_index[2],colour_index[1],colour_index[0]);
              fprintf(fp,"}}\n");
              

/* CURRENT - attempting to fix up Josh's dodgy colour stuff :-) */
 /* 
              switch (style)
                {
                case HSV:
                  fprintf(fp,"texture_list{HSV%d HSV%d HSV%d}}\n",
                              colour_index[2],colour_index[1],colour_index[0]);
                  break;

                case REDWHITEBLUE:
                  fprintf(fp,"texture_list { ");
                  fprintf(fp,"RGB_%d RGB_%d RGB_%d ",
                              colour_index[2],colour_index[1],colour_index[0]);
                  fprintf(fp,"}}\n");
                  break;
                }
  */
              }
            fprintf(fp,"}\nobject { Surface translate <0,0,0> rotate <0,0,0> }\n");
            break;

          case GL_TRIANGLE_STRIP:

            p1 = rlist->data; 
            ARR3SET(v1, vec2);
            ARR3ADD(v1, p1->rx);
            ARR3SET(n1, p1->rn);

            rlist = g_slist_next(rlist);

            p2 = rlist->data;
            ARR3SET(v2, vec2);
            ARR3ADD(v2, p2->rx);
            ARR3SET(n2, p2->rn);

            rlist = g_slist_next(rlist);

            while (rlist)
              {
              p3 = rlist->data;
              ARR3SET(v3, vec2);
              ARR3ADD(v3, p3->rx);
              ARR3SET(n3, p3->rn);

              fprintf(fp, "triangle { \n");

              fprintf(fp, "< %f, %f, %f>", v1[0], v1[1], v1[2]);
              fprintf(fp, "< %f, %f, %f>", v2[0], v2[1], v2[2]);
              fprintf(fp, "< %f, %f, %f>\n", v3[0], v3[1], v3[2]);
 
/* FIXME - can we have individual colours at each vertex? */
              fprintf(fp, " texture { \n");
              fprintf(fp, " pigment { color rgbt<%f, %f, %f, %f> } }\n",
                           p3->colour[0], p3->colour[1], p3->colour[2],
                           1.0-sysenv.render.transmit);

/* NB: POVRay highlight size has reverse sense. */
			  fprintf(fp, " finish{Phong_Shiny}}}\n");
			  /*
              fprintf(fp, " finish { phong %f phong_size %d } }\n",
                           sysenv.render.shl_strength, (gint) sysenv.render.shl_size);*/

              ARR3SET(v1, v2);
              ARR3SET(v2, v3);
              ARR3SET(n1, n2);
              ARR3SET(n2, n3);

              rlist = g_slist_next(rlist);
              }
            break;

/* general case */
          default:
            fprintf(fp, "polygon { %d\n", g_slist_length(rlist));
            while (rlist)
              {
              p1 = rlist->data;
              ARR3SET(vec1, vec2);
              ARR3ADD(vec1, p1->rx);
              fprintf(fp, "< %f, %f, %f>", vec1[0], vec1[1], vec1[2]);
  
              rlist = g_slist_next(rlist);
  
              if (rlist)
                fprintf(fp, ",\n");
              else
                fprintf(fp, "\n");
              }

		/* FIXME - can we have individual colours at each vertex? */
        fprintf(fp, " texture { \n");
        fprintf(fp, " pigment { color rgbt<%f, %f, %f, %f> }\n",
                                               p1->colour[0],
                                               p1->colour[1],
                                               p1->colour[2],
                                      1.0-sysenv.render.transmit);
			/* NB: POVRay highlight size has reverse sense. */
			fprintf(fp, " finish{Phong_Shiny}}}\n");
		/*
        fprintf(fp, " finish { phong %f phong_size %d } }\n",
          sysenv.render.shl_strength, (gint) sysenv.render.shl_size); */

          }
        }
      while (plist);
      break;
    }
  }

fclose(fp);
return(0);
}

