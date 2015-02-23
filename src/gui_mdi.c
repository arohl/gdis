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
#include <stdlib.h>
#include <string.h>

#ifndef __WIN32
#include <sys/times.h>
#endif

#include "gdis.h"
#include "file.h"
#include "coords.h"
#include "edit.h"
#include "zone.h"
#include "matrix.h"
#include "model.h"
#include "render.h"
#include "select.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "mdi_pak.h"

#define DEBUG 0

/****************/
extern struct sysenv_pak sysenv;
/****************/

void mdi_model_create(gint);
void rot(gdouble *, gdouble , gdouble , gdouble);
GtkWidget *spinner[MAX_MODELS];
struct mdi_pak mdi_data;

/*******************/
/* MDI setup & run */
/*******************/
void mdi_model_setup(void)
{
GtkSpinButton *spin;
gint i,num_req,last_spin;
gdouble biggest;
GSList *list;
struct model_pak *data;

/* box size required */
spin = GTK_SPIN_BUTTON (spinner[0]);
mdi_data.box_dim = gtk_spin_button_get_value_as_int(spin);

/* NB: if an MDI run has already been done - terminate on this model */
last_spin = 0;
for (list=sysenv.mal ; list ; list=g_slist_next(list))
  {
  data = list->data;
  if (data->id == MDI)
    break;
  last_spin++;
  }

/* numer of loaded models = max possible components */
mdi_data.num_comp = last_spin;
mdi_data.comp_idx = g_malloc0(mdi_data.num_comp*sizeof(gint));
mdi_data.comp_req = g_malloc0(mdi_data.num_comp*sizeof(gint));
mdi_data.comp_done = g_malloc0(mdi_data.num_comp*sizeof(gint));

/* error trap - don't use an MDI model as a (non-zero) component */
if (!last_spin)
  {
  gui_text_show(ERROR, "Have you done something silly?\n");
  return;
  }

/* go through non solvent molecules & sum amount required */
num_req=0;
for (i=1 ; i<last_spin ; i++)
  {
  spin = GTK_SPIN_BUTTON(spinner[i]);
  mdi_data.comp_req[i] = gtk_spin_button_get_value_as_int(spin);
  num_req += mdi_data.comp_req[i];
/* redundant? */
  mdi_data.comp_idx[i] = i;
  }
/* calculate for solvent */
mdi_data.comp_req[0] = pow(mdi_data.box_dim,3) - num_req;
/* 8 corner sites - must have solvent */
/* TODO - is this really necessary??? */
if (mdi_data.comp_req[0] < 9)
  {
  gui_text_show(ERROR, "Too many solute molecules required!\n");
  return;
  }
mdi_data.comp_idx[0] = 0;
mdi_data.comp_done[0] = 0;

/* auto calculate lattice spacing */
biggest=-1.0;
for(i=0 ; i<mdi_data.num_comp ; i++)
  if (mdi_data.comp_req[i])
    {
    data = model_ptr(i, RECALL); 
    g_return_if_fail(data != NULL);
    if (data->rmax > biggest)
      biggest = data->rmax;
    }

/* auto setting = 2*biggest model's radius + safety */
mdi_data.latt_sep=2.0*biggest + 1.0; 

#if DEBUG
printf("Boxdim = %d\n",mdi_data.box_dim);
printf("Lattice sep = %f\n",mdi_data.latt_sep);
for (i=0 ; i<mdi_data.num_comp ; i++)
  {
  printf("Model %d, idx=%d, req=%d, done=%d\n",i,mdi_data.comp_idx[i]
                                                ,mdi_data.comp_req[i]
                                                ,mdi_data.comp_done[i]);
  }
#endif

/* main call */
mdi_model_create(last_spin);

/* cleanup */
g_free(mdi_data.comp_idx);
g_free(mdi_data.comp_req);
g_free(mdi_data.comp_done);
}

/****************************/
/* GROMACS genbox solvation */
/****************************/
/* return true on succes */
gint gromacs_genbox(void)
{
gchar *genbox, *cmd, *inp, *out;
struct model_pak *model;

model = sysenv.active_model;
if (!model)
  return(FALSE);

genbox = file_find_program("genbox");
if (!genbox)
  return(FALSE);

printf("using: %s\n", genbox);

inp = gun("gro");
out = gun("gro");

write_gromacs(inp, model);

cmd = g_strdup_printf("%s -cs -cp %s -o %s", genbox, inp, out);

system(cmd);

/* TODO - delete old atoms */
core_delete_all(model);

read_gromacs_gro(out, model);

zone_init(model);

core_render_mode_set(STICK, model->cores);

g_free(cmd);

g_free(inp);
g_free(out);

g_free(genbox);

return(TRUE);
}

/********************/
/* MDI input dialog */
/********************/
void gui_mdi_dialog(void)
{
gint i;
gpointer dialog;
GtkWidget *window, *vbox, *label, *frame;
GtkAdjustment *adj;
GString *frame_label;
GSList *list;
struct model_pak *data;

/* NEW - hijack if genbox (GROMACS solvator) is found */
if (gromacs_genbox())
  return;

/* we need some models to work with */
if (g_slist_length(sysenv.mal) < 2)
  {
  gui_text_show(ERROR, "Please load solvent and solute molecules first.\n"); 
  return;
  }

/* get solvent */
data = model_ptr(0, RECALL);
g_return_if_fail(data != NULL);

/* try to get a dialog id, associate with solvent model */
dialog = dialog_request(MDI, "MD initializer", NULL, NULL, data);
if (!dialog)
  return;
window = dialog_window(dialog);

/* create input (spin) widgets */
/* box dimensions */
frame_label = g_string_new(NULL);
g_string_sprintf(frame_label,"Box dimension");
frame = gtk_frame_new(frame_label->str);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, TRUE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
 
label = gtk_label_new ("Side length (lattice points)");
/*  gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5); */
gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, TRUE, 0);

adj = (GtkAdjustment *) gtk_adjustment_new (4, 3, 100, 1, 1, 0);

spinner[0] = gtk_spin_button_new (adj, 0, 0);
gtk_spin_button_set_wrap(GTK_SPIN_BUTTON (spinner[0]), FALSE);
gtk_box_pack_start(GTK_BOX (vbox), spinner[0], FALSE, TRUE, 0);

g_string_sprintf(frame_label,"Model: %s",data->basename);

frame = gtk_frame_new(frame_label->str);
gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, TRUE, 0);
vbox = gtk_vbox_new(FALSE, 0);
gtk_container_set_border_width(GTK_CONTAINER (vbox), PANEL_SPACING);
gtk_container_add(GTK_CONTAINER(frame), vbox);
label = gtk_label_new ("This will be the solvent");
gtk_box_pack_start(GTK_BOX(vbox), label, FALSE, TRUE, 0);

/* omit solvent (model 0) */
i = 1;
list = g_slist_next(sysenv.mal);
while (list)
  {
  data = list->data;

/* don't include any previous MDI model */
  if (data->id == MDI)
    continue;
  g_string_sprintf(frame_label,"Model: %s",data->basename);

  frame = gtk_frame_new(frame_label->str);
  gtk_box_pack_start(GTK_BOX(GTK_DIALOG(window)->vbox), frame, TRUE, TRUE, 0);
  vbox = gtk_vbox_new(FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(vbox), PANEL_SPACING);
  gtk_container_add(GTK_CONTAINER(frame), vbox);
 
  label = gtk_label_new ("Number required");
/*  gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5); */
  gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, TRUE, 0);

  adj = (GtkAdjustment *) gtk_adjustment_new (0, 0, 100, 1, 5, 0);

  spinner[i] = gtk_spin_button_new (adj, 0, 0);
  gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (spinner[i]), TRUE);
  gtk_box_pack_start (GTK_BOX (vbox), spinner[i], FALSE, TRUE, 0);

  list = g_slist_next(list);
  i++;
  }

/* buttons */
gui_stock_button(GTK_STOCK_EXECUTE, mdi_model_setup, NULL,
                   GTK_DIALOG(window)->action_area);

gui_stock_button(GTK_STOCK_CLOSE, dialog_destroy, dialog,
                   GTK_DIALOG(window)->action_area);


/* display the dialog */
gtk_widget_show_all(window);

g_string_free(frame_label, TRUE);
}

/***********/
/* run MDI */
/***********/
#define DEBUG_MDI_MODEL_CREATE 0
void mdi_model_create(gint new)
{
gint i,j,na,nb,ns,replace,ret;
struct model_pak *data;
#ifndef __WIN32
struct tms buffer;
#endif

/* main array */
mdi_data.array = (gint *) g_malloc0(mdi_data.box_dim*mdi_data.box_dim*
                                   mdi_data.box_dim*sizeof(gint));

/* calculate the total number of atoms & bonds & sites in the box */
na = nb = 0;
ns = -mdi_data.comp_req[0];
for (i=0 ; i<mdi_data.num_comp ; i++)
  {
  data = model_ptr(i, RECALL);
  mdi_data.comp_done[i]=0;
  }

replace=1;

data = model_new();

/* check & init */
g_return_if_fail(data != NULL);
data->id = MDI;
strcpy(data->filename,"MDI_model");
g_free(data->basename);
data->basename = g_strdup("MDI model");


/* initialize random generator */
#if __WIN32
j=666;
#else
j=times(&buffer);
#endif

#if DEBUG_MDI_MODEL_CREATE 
printf("Random seed = %d\n",j);
#endif
srand(j);

/* do the business */
ret = fill();

#if DEBUG_MDI_MODEL_CREATE 
printf("fill() return code: %d\n",ret);
#endif
if (ret == 2)
  {
  gui_text_show(ERROR, "Run was unsuccessful!\n");

  model_delete(data);
  g_free(mdi_data.array);
  return;
  }

/* copy box into the model coord data array */
write_dat(data);
model_prep(data);

/* update model tree display */
if (replace)
  {
  tree_model_add(data);
  gui_model_select(data);
  }
else
  {
  tree_model_add(data);
  }

#if DEBUG_MDI_MODEL_CREATE 
printf("done model init\n");
#endif

g_free(mdi_data.array);

#if DEBUG_MDI_MODEL_CREATE 
printf("MDI done.\n");
#endif
}

/*********************/
/* Main box fill sub */
/*********************/
gint fill()
{
gint num_sites, max_sites,num_cand;
gint i,j,k,di,dj,dk,r2,s;
gint req_tot,ci,ri,pos;
gint *r2_min, max_min;
gdouble ran2();
struct cand_pak *cand;
struct box_pak *site;

/* allocate for sites - ie solvent not included XCPT corners */
max_sites = pow(mdi_data.box_dim,3) - mdi_data.comp_req[0] + 8;
site = g_malloc(max_sites*sizeof(struct box_pak));

r2_min = (gint *) g_malloc(pow(mdi_data.box_dim,3)*sizeof(gint));

/* max out the number of candidate sites */
cand = g_malloc(pow(mdi_data.box_dim,3) *sizeof(struct cand_pak));

/* NB: initially every point in box should be 0 => solvent */
#if DEBUG
printf("Filling box: %d x %d x %d\n",mdi_data.box_dim,
                                     mdi_data.box_dim,mdi_data.box_dim);
#endif

/*  NO cadidate inititalization */
num_sites = 0;

/* compute total number of components to dissolve */
req_tot=0;
for (i=1 ; i<mdi_data.num_comp ; i++)    /* omit solvent */
  req_tot += mdi_data.comp_req[i];

ci = 1;                         /* starting component number */

/* loop over total number to dissolve */
for (ri=0 ; ri<req_tot ; ri++)
  {
/* seek a component which doesn't yet have the required number solvated */
  while(ci < g_slist_length(sysenv.mal) && 
        mdi_data.comp_done[ci] >= mdi_data.comp_req[ci])
    ci++; 
  if (ci == g_slist_length(sysenv.mal))
    {
    printf("fill() error: mismatch in component requirements!\n");
    return(2);
    }
  mdi_data.comp_done[ci]++;

/* evaluate all sites - sphere of clearance */
pos=0;
max_min=0;
for(i=0 ; i<mdi_data.box_dim ; i++)
  {
  for(j=0 ; j<mdi_data.box_dim ; j++)
    {
    for(k=0 ; k<mdi_data.box_dim ; k++)
      {
/* find minimum distance from this site to all candidate sites */
      r2_min[pos] = pow(mdi_data.box_dim,3);
      for (s=0 ; s<num_sites ; s++)
        {
if (s == max_sites)
  {
  printf("Error: site indexing beyond boundary!\n");
  return(2);
  }
/* PBC */
        di = abs(site[s].x-i);
        if (di > mdi_data.box_dim/2.0) di -= mdi_data.box_dim/2.0;
        dj = abs(site[s].y-j);
        if (dj > mdi_data.box_dim/2.0) dj -= mdi_data.box_dim/2.0;
        dk = abs(site[s].z-k);
        if (dk > mdi_data.box_dim/2.0) dk -= mdi_data.box_dim/2.0;

        r2 = di*di + dj*dj + dk*dk;
/* test for new minimum */
        if (r2 < r2_min[pos])
          r2_min[pos] = r2;
        }
/* get the overall maximum (of the minima) */
      if (r2_min[pos] > max_min) max_min = r2_min[pos];
      pos++;
      }
    }
  }

s=pos=0;
for(i=0 ; i<mdi_data.box_dim ; i++)
  {
  for(j=0 ; j<mdi_data.box_dim ; j++)
    {
    for(k=0 ; k<mdi_data.box_dim ; k++)
      {
      if (r2_min[pos] == max_min) 
        {
        cand[s].pos = pos;
        cand[s].x = i; 
        cand[s].y = j; 
        cand[s].z = k; 
        s++;
        }
      pos++;
      }
    }
  }

num_cand = s;
if (!num_cand)
  {
  printf("No candidate sites found!\n");
  g_free(r2_min);
  g_free(site);
  g_free(cand);
  return(2);
  }
#if DEBUG
printf("Found %d candidate site(s)\n",num_cand);
#endif

/* Select (at random) a candidate site into */
/* which the current component will be placed */

  s = (gint) (ran2() * num_cand);

#if DEBUG
printf("s=%d\n",s);
#endif

/* do site update */
  site[num_sites].component = ci;
  site[num_sites].x = cand[s].x;
  site[num_sites].y = cand[s].y;
  site[num_sites].z = cand[s].z;
  num_sites++;

/* do array data update */
  pos = cand[s].pos;
  *(mdi_data.array+pos) = ci; 
  }

/* final data in nice array format */
#if DEBUG
pos=0;
for (i=0; i<mdi_data.box_dim ; i++)
  {
  for (j=0; j<mdi_data.box_dim ; j++)
    {
    for (k=0; k<mdi_data.box_dim ; k++)
      {
      printf("%d",*(mdi_data.array+pos));
      pos++;
      }
    printf("\n");
    }
  printf("\n");
  }
printf("fill() done\n");
#endif

g_free(r2_min);
g_free(site);
g_free(cand);
return(0);
}

/*********************/
/* Save box to model */
/*********************/
void write_dat(struct model_pak *dest)
{
gint s,m,pos,dat;
gdouble x,y,z,phi,theta,psi,start,stop,step;
gdouble q[3],box_len;
gdouble ran2();
GSList *list;
struct model_pak *data;
struct core_pak *core, *copy;

#if DEBUG
printf("Starting write_dat()...\n");
#endif

box_len = mdi_data.box_dim * mdi_data.latt_sep;

start = mdi_data.latt_sep/2.0;
stop = box_len;
step = mdi_data.latt_sep;

#if DEBUG
printf("Box dim = %f\n",box_len);
printf("Box rng = %f,%f,%f\n",start,stop,step);
#endif

pos=0; /* pointer - current array position */
dat=0; /* new model - atom coords pointer */


for (z=start ; z<stop ; z+=step)
  {
  for (y=start ; y<stop ; y+=step)
    {
    for (x=start ; x<stop ; x+=step)
      {
/* trap */
if (pos >= pow(mdi_data.box_dim,3))
  {
  printf("Error at (%f,%f,%f)\n",x,y,z);
  printf("Program bug: address (pos) out of bounds!\n");
  return;
  }

/* get component at this (x,y,z) pos'n */
      s = *(mdi_data.array+pos);
/* get the model this corresponds to */
      m = mdi_data.comp_idx[s]; 
      data = model_ptr(m, RECALL);
      pos++;
/* create some random rotations */
      phi = 2.0*PI*ran2();
      theta = PI - 2.0*PI*ran2();
      psi = 2.0*PI*ran2();

      for (list=data->cores ; list ; list=g_slist_next(list))
        {
        core = list->data;
        if (core->status & DELETED)
          continue;
 
/* do rotation */
        ARR3SET(q, core->rx);
        rot(&q[0],phi,theta,psi);
/* do lattice site translation of coords */
        q[0] += x;
        q[1] += y;
        q[2] += z;

/* write the atom data to the model structure */
        copy = dup_core(core);
        dest->cores = g_slist_prepend(dest->cores, copy);

        copy->status = NORMAL;
        copy->orig = TRUE;
        copy->primary = TRUE;
        copy->primary_core = NULL;
        copy->region = REGION1A;
        ARR3SET(copy->x, q);
        }
      }
    }
  }

/* model info */
dest->id = MDI;
dest->fractional = FALSE;
dest->periodic = 3;
dest->pbc[0] = box_len;
dest->pbc[1] = box_len;
dest->pbc[2] = box_len;
dest->pbc[3] = PI/2.0;
dest->pbc[4] = PI/2.0;
dest->pbc[5] = PI/2.0;

#if DEBUG
printf("write_dat() done\n");
#endif
}

/***********************************************/
/* input a 3 element array = coords to rotate  */
/*       plus angles to rotate by              */
/* output the three coords rotated about 0,0,0 */
/*        by input angles                      */
/***********************************************/
void rot(gdouble *coords, gdouble phi, gdouble theta, gdouble psi)
{
gdouble sph, cph, sth, cth, sps, cps;
gdouble rot[3][3];

/* let's make things easier */
sph=sin(phi);
cph=cos(phi);
sth=sin(theta);
cth=cos(theta);
sps=sin(psi);
cps=cos(psi);

/* make the matrix */
rot[0][0]= cps*cph-cth*sph*sps;
rot[0][1]= cps*sph+cth*cph*sps;
rot[0][2]= sps*sth;
rot[1][0]=-sps*cph-cth*sph*cps;
rot[1][1]=-sps*sph+cth*cph*cps;
rot[1][2]= cps*sth;
rot[2][0]= sth*sph;
rot[2][1]=-sth*cph;
rot[2][2]= cth;

/* multiply the coords with the matrix */
vecmat(&rot[0][0],coords);
}

/*********************/
/* My random routine */
/*********************/
gdouble ran2()
{
gint i,j;
gdouble f;

i=j=0;
/* generate 2 random (non 0) integers */
while(!i)
  i = rand();
while(!j)
  j = rand();
/* convert to a single float (0.0,1.0] */
if (i > j)
  f = (gdouble) j / (gdouble) i;
else
  f = (gdouble) i / (gdouble) j;

return f;
}

