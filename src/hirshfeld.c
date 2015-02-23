/*
 *  hirshfeld.c
 *  newsurface
 *
 *  Created by Anthony S. Mitchell and Joshua J. McKinnon.
 *  Contact:  Joshua.McKinnon@une.edu.au
 *  Copyright (c) 2004 University of New England. All rights reserved.
 *
 * This file contains (almost) all of the routines required to calculate
 * w(r) for the Hirshfeld surface.
 *
 
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
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <signal.h>
#include <time.h>

#include "gdis.h"
#include "coords.h"
#include "matrix.h"
#include "space.h"
#include "numeric.h"
#include "hirshfeld.h"
#include "hirshfeld_data.h"
#include "molsurf.h"
#define MAXATOMS 54
#define NPOINTS 700
/* spline limits (in Angstroms) */
#define MAXR 7.4084815
/* crystal environment cutoff */
//#define HFS_RMAX 10.0
/* JJM DEBUG try a smaller cutoff */
#define HFS_RMAX 5.0

gboolean hfs_rho_init = FALSE;
struct rho_pak
{
gdouble x[NPOINTS];
gdouble y[NPOINTS];
gdouble y2[NPOINTS];
} atomdata[MAXATOMS];

/************************************************************/
/* initialization routine for Hirshfeld surface calculation */
/************************************************************/
void hfs_init(void)
{
gint i, j;
gdouble x /* , delta, big, small */ ;

/* only compute spline once */
if (hfs_rho_init)
  return;

/* setup data array for spline fit */
for (i=0 ; i<MAXATOMS ; i++)
  {
  for (j=0 ; j<NPOINTS ; j++)
    {
    x = (gdouble) j * MAXR / (gdouble) NPOINTS;		
/*** JJM: evaldens at sqrt(r) ***/
    atomdata[i].x[j] = x*x;
    atomdata[i].y[j] = rho[i][j];
    }

/***  use a crude numerical derivative at the first point  ***/
/*
  delta = 0.0001;
  big = evaldens(atom[i],0.0) - evaldens(atom[i],sqrt(delta));
  big = big/delta;
*/
 
/***  use a crude numerical derivative at the last point  ***/
/*
  delta = 0.0001;
  small = evaldens(atom[i],sqrt(maxr)) - evaldens(atom[i],(sqrt(maxr)+delta));
  small = small/delta;
*/

/* FIXME - I've just used 0.0 for the end point derivatives */
  spline(atomdata[i].x, atomdata[i].y, NPOINTS, 0.0, 0.0, atomdata[i].y2);
  }

hfs_rho_init = TRUE;
}

/*******************************************************************/
/* return a list of cores (ie molecule) closest to the given point */
/*******************************************************************/
GSList *hfs_closest_molecule(gdouble *x, struct model_pak *model)
{
gdouble dr, min, r[3];
GSList *list;
struct core_pak *core;
struct mol_pak *mol=NULL;

/*
printf("closest to");
P3VEC(" : ", x);
*/

g_assert(model != NULL);

min = 1000000000.0;

/*
printf("min = %f\n", min);
*/

for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  ARR3SET(r, core->x);
  vecmat(model->latmat, r);

/*
P3VEC("core: ", r);
*/

  ARR3SUB(r, x);

  dr = VEC3MAGSQ(r);
  if (dr < min)
    {
    mol = core->mol;
    min = dr;
    }
  }
/*
printf("min = %f\n", min);
*/

if (mol)
  return(mol->cores);
return(NULL);
}

/***********************************/ 
/* generate atoms around selection */
/***********************************/ 

GSList *hfs_calc_env(GSList *select, struct model_pak *model)
{
gboolean flag, primary;
gint i, imax[3];
gdouble r2, dx[3];
GSList *ilist, *list, *list2, *nlist;
struct image_pak *image;
struct core_pak *core, *copy, *test;

/* checks */
g_assert(model != NULL);

/* construct periodic image vectors */
for (i=0 ; i<model->periodic ; i++)
  imax[i] = HFS_RMAX/model->pbc[i] + 1;

model->image_limit[0] = imax[0];
model->image_limit[1] = imax[0]+1;
model->image_limit[2] = imax[1];
model->image_limit[3] = imax[1]+1;
model->image_limit[4] = imax[2];
model->image_limit[5] = imax[2]+1;
space_make_images(CREATE, model);
coords_compute(model);

/* original + image iteration */
nlist = NULL;
for (list=model->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* mark cores in the central molecule */
  if (g_slist_find(select, core))
    primary = TRUE;
  else
    primary = FALSE;

/* add the atom */
  ilist = NULL;
  do
    {
/* add coordinates */
/* NB: probably don't need to worry about shells, as this is */
/* a temporary model and the algorithm will ignore them anyway */
    copy = dup_core(core);
    if (ilist)
      {
/* image translation */
      image = ilist->data;
      ARR3ADD(copy->rx, image->rx);
      ilist = g_slist_next(ilist);
      }
    else
      ilist = model->images;

/* set cartesian coords */
    ARR3SET(copy->x, copy->rx);
/* JJM added line below so that copy->x contains fractional coords.
    Should I take account of the centroid here? */
    vecmat(model->ilatmat,copy->x);
    ARR3ADD(copy->x, model->centroid); 
    
    copy->primary = primary;

/* test if core satisfies distance condition */
    flag = FALSE;
    for (list2=select ; list2 ; list2=g_slist_next(list2))
      {
      test = list2->data;

      ARR3SET(dx, test->rx);
      ARR3SUB(dx, copy->rx);
      r2 = VEC3MAGSQ(dx);
      if (r2 < HFS_RMAX*HFS_RMAX)
        {
        flag = TRUE;
        break;
        }
      }
    if (flag)
      nlist = g_slist_prepend(nlist, copy);
    }
  while (ilist);
  }

/* reset source model's images */
space_make_images(INITIAL, model);

#if DEBUG_HFS_CALC_ENV
printf("hfs_calc_env: %d cores\n", g_slist_length(nlist));
#endif

return(nlist);
}

/************************************/
/* retrieve an interpolated density */
/************************************/
gdouble hfs_calc_density(gint n, gdouble r2)
{
gdouble value;


/* FIXME - no support for atoms up to (and not including) Ba */
/*
g_assert(n < MAXATOMS);
*/
n = CLAMP(n, 1, MAXATOMS);

if (r2 > atomdata[n].x[NPOINTS-1])
  return(0.0);

splint(atomdata[n].x, atomdata[n].y, atomdata[n].y2, NPOINTS, r2, &value);
return(value);
}

/********************************************/
/* Use the selection to calculate        	*/
/* sensible (but not foolproof!!) limits 	*/
/* for the marching cubes search, of the	*/
/* extremeties of the selection atoms in	*/
/* x,y,z plus <tolerance> in each direction */
/********************************************/
int hfs_calulation_limits(GSList *selection, struct model_pak *model, gdouble *min, gdouble *max)
{
    
    gdouble tolerance[3] = {3.0, 3.0, 3.0}; /* Angstrom */
    GSList *list;
    struct core_pak *thisAtom;
    int i;
    for (i=0; i<3; i++){
        min[i] = 3e10;
        max[i] =-3e10;
    }
    vecmat(model->ilatmat,tolerance);
    for (list = selection; list; list=g_slist_next(list)){
        thisAtom = list->data;
        for (i=0; i<3; i++){
          if ( thisAtom->x[i] < min[i]) min[i] = thisAtom->x[i];
          if ( thisAtom->x[i] > max[i]) max[i] = thisAtom->x[i];   
        }
    }
    for (i=0; i<3; i++) {
        min[i] = min[i] - tolerance[i];
        max[i] = max[i] + tolerance[i];
    }
    return(0);
}

/*************************************/
/* Use the gradient of the Hirshfeld */
/* surface to calculate the normal   */
/*************************************/

int hfs_calc_normals(GSList *points, struct model_pak *model, GSList *myThing)
{
    gdouble delta = 5.0e-3;
    gdouble gx, gy, gz; 
    gdouble wPoint, wgx, wgy, wgz;
    gdouble xDelta[3],yDelta[3],zDelta[3];
    gdouble temp[3];
    GSList *pointList;
    struct smv_pak *thisPoint;
    time_t jobStartTime, jobEndTime;
    
    jobStartTime = time(NULL);
    delta = 5.0e-3;
    for (pointList=points; pointList; pointList=g_slist_next(pointList)){
        
        thisPoint = pointList->data;
        if (!thisPoint->adj)
            continue;

        ARR3SET(temp, thisPoint->rx);
               
        gx = thisPoint->rx[0] + delta;
        gy = thisPoint->rx[1] + delta;
        gz = thisPoint->rx[2] + delta;
        VEC3SET(xDelta,gx,thisPoint->rx[1],thisPoint->rx[2]);
        VEC3SET(yDelta,thisPoint->rx[0],gy,thisPoint->rx[2]);
        VEC3SET(zDelta,thisPoint->rx[0],thisPoint->rx[1],gz);
        wPoint  = hfs_calc_wf(thisPoint->rx[0],thisPoint->rx[1],thisPoint->rx[2],model,myThing);
        wgx = hfs_calc_wf(xDelta[0],xDelta[1],xDelta[2],model,myThing);
        wgy = hfs_calc_wf(yDelta[0],yDelta[1],yDelta[2],model,myThing);
        wgz = hfs_calc_wf(zDelta[0],zDelta[1],zDelta[2],model,myThing);
        
        thisPoint->nx[0] = (wgx - wPoint)/delta;
        thisPoint->nx[1] = (wgy - wPoint)/delta;
        thisPoint->nx[2] = (wgz - wPoint)/delta;
//        normalize(thisPoint->nx,3);
        ARR3SET(thisPoint->n, thisPoint->nx);
//        P3VEC("Normal:",thisPoint->nx);
        
   }
   jobEndTime = time(NULL);
#if DEBUG_HFS_NORMALS
   fprintf(stderr," [%.1fs]\n",(float)(jobEndTime - jobStartTime));
#endif
   return(0);
}

/*************************************/
/* Use the gradient of the promolecule */
/* surface to calculate the normal   */
/*************************************/

int ssatoms_calc_normals(GSList *points, struct model_pak *model)
{
    gdouble delta = 5.0e-3;
    gdouble gx, gy, gz; 
    gdouble wPoint, wgx, wgy, wgz;
    gdouble xDelta[3],yDelta[3],zDelta[3];
    gdouble temp[3];
    GSList *pointList;
    struct smv_pak *thisPoint;
    time_t jobStartTime, jobEndTime;
    
    jobStartTime = time(NULL);
    delta = 5.0e-3;
    for (pointList=points; pointList; pointList=g_slist_next(pointList)){
        
        thisPoint = pointList->data;
        if (!thisPoint->adj)
            continue;
		
        ARR3SET(temp, thisPoint->rx);
		
        gx = thisPoint->rx[0] + delta;
        gy = thisPoint->rx[1] + delta;
        gz = thisPoint->rx[2] + delta;
        VEC3SET(xDelta,gx,thisPoint->rx[1],thisPoint->rx[2]);
        VEC3SET(yDelta,thisPoint->rx[0],gy,thisPoint->rx[2]);
        VEC3SET(zDelta,thisPoint->rx[0],thisPoint->rx[1],gz);
        wPoint  = ssatoms_calc_density(thisPoint->rx[0],thisPoint->rx[1],thisPoint->rx[2],model);
        wgx = ssatoms_calc_density(xDelta[0],xDelta[1],xDelta[2],model);
        wgy = ssatoms_calc_density(yDelta[0],yDelta[1],yDelta[2],model);
        wgz = ssatoms_calc_density(zDelta[0],zDelta[1],zDelta[2],model);
        
        thisPoint->nx[0] = (wgx - wPoint)/delta;
        thisPoint->nx[1] = (wgy - wPoint)/delta;
        thisPoint->nx[2] = (wgz - wPoint)/delta;
		normalize(thisPoint->nx,3);
        ARR3SET(thisPoint->n, thisPoint->nx);
//		P3VEC("Normal:",thisPoint->nx);
        
	}
	jobEndTime = time(NULL);
#if DEBUG_HFS_NORMALS
	fprintf(stderr," [%.1fs]\n",(float)(jobEndTime - jobStartTime));
#endif
	return(0);
}

/******************************************/
/* compute the weight function at a point */
/*					  */
/* Note we now pass three doubles (x,y,z) */
/* instead of a single vector because it  */
/* makes life much easier for the second  */
/* derivative calculations 		  */
/******************************************/
gdouble hfs_calc_wf(gdouble x, gdouble y, gdouble z, struct model_pak *model, GSList *myThing)
{
gdouble r2;
gdouble weight_r=0.0, r[3];
gdouble point[3]; /*vector for the point (x,y,z) */
gdouble rhoNugget, rhoMol, rhoXtal;
GSList *list, *mlist;
struct core_pak *core;

VEC3SET(point,x,y,z);
/* checks */
g_assert(model != NULL);
g_assert(hfs_rho_init == TRUE);
/* add up the density contribution from the molecule */

/*** JJM The 'molecule' that we're calculating the Hirshfeld surface of should be the selected molecule. Hence the new definition of mlist here ***/

mlist = model->selection;



rhoMol=0.0;

/***J loop over the atoms in the selected molecule J***/
for (list=mlist ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* compute cartesian distance squared (in Angs) */
//  vecmat(model->latmat,core->x);
  ARR3SET(r, core->rx);
  vecmat(model->ilatmat,r);
  ARR3ADD(r, model->centroid);
  vecmat(model->latmat,r);
  ARR3SUB(r, point);
  r2 = VEC3MAGSQ(r);

/*  ARR3SET(coreCartesian, core->x);
  vecmat(model->latmat,coreCartesian);
  ARR3SUB(coreCartesian,x);
  r2 = VEC3MAGSQ(r);
*/  
  rhoNugget = hfs_calc_density(core->atom_code, r2);
  rhoMol += rhoNugget;
  }
/* add up the total density */
rhoXtal=0.0;

/***JJM loop over the atoms in the crystal J***/
for (list=myThing ; list ; list=g_slist_next(list))
  {
  core = list->data;

/* compute cartesian distance squared (in Angs) */
//  vecmat(model->latmat,core->x);
  ARR3SET(r, core->rx);
    vecmat(model->ilatmat,r);
    ARR3ADD(r, model->centroid);
    vecmat(model->latmat,r);
    
#ifdef JJM_DEBUG_EXTREME02
P3VEC(" r: ", r);
P3VEC(" x: ", point);
#endif
  ARR3SUB(r, point);
  r2 = VEC3MAGSQ(r);

/*  ARR3SET(coreCartesian, core->x);
  vecmat(model->latmat,coreCartesian);
  ARR3SUB(coreCartesian,x);
  r2 = VEC3MAGSQ(r);
  */  
  rhoNugget = hfs_calc_density(core->atom_code, r2);
  rhoXtal += rhoNugget;
  }
#ifdef JJM_DEBUG02
if ( rhoMol > 0.0 || rhoXtal > 0.0){
    fprintf(stderr,"rhoMol: %.3f ___ rhoXtal: %.3f\n",rhoMol,rhoXtal);
}
#endif

/* calc ratio */
if (rhoXtal != 0.0){

/* NEW */
  if (rhoXtal > rhoMol)
    weight_r = rhoMol/rhoXtal;
    else {
        fprintf(stderr,"WARNING: Weight function > 1.0. Value: %f\n",weight_r);
 	weight_r = 0.0;
    }
    
}
else
  weight_r = 0.0;

#define STRICT 1
#if STRICT
g_assert(weight_r <= 1.0);
#else
if (weight_r > 1.0)
  weight_r = 1.0;
#endif
/* JJM debug
if (weight_r > 0.5){
    fprintf(stderr,"w %.2f  ",weight_r);
    P3VEC("at ",x);
    } 
*/
return(weight_r);
}

/*********************************************************************/
/* compute the sum of spherical atom density at a point 			 */
/* a.k.a. Promolecule												 */
/* See Mitchell & Spackman, J. Comp. Chem., v. 21, pp 933-942 (2000) */
/*																	 */
/*																	 */
/*					  												 */
/* Note we now pass three doubles (x,y,z) 							 */
/* instead of a single vector because it  							 */
/* makes life much easier for the second  							 */
/* derivative calculations 		  									 */
/*********************************************************************/
gdouble ssatoms_calc_density(gdouble x, gdouble y, gdouble z, struct model_pak *model)
{
	gdouble r2;
	gdouble r[3];
	gdouble point[3]; /*vector for the point (x,y,z) */
	gdouble rhoNugget, rhoMol;
	GSList *list, *mlist;
	struct core_pak *core;
	
	VEC3SET(point,x,y,z);
	/* checks */
	g_assert(model != NULL);
	g_assert(hfs_rho_init == TRUE);
	/* add up the density contribution from the molecule */
	
	/*** JJM The 'molecule' that we're calculating the  surface of should be the selected molecule. Hence the new definition of mlist here ***/
	
	mlist = model->selection;
	
	
	
	rhoMol=0.0;
	
	/***J loop over the atoms in the selected molecule J***/
	for (list=mlist ; list ; list=g_slist_next(list))
	{
		core = list->data;
		
		/* compute cartesian distance squared (in Angs) */
		//  vecmat(model->latmat,core->x);
		ARR3SET(r, core->rx);
		vecmat(model->ilatmat,r);
		ARR3ADD(r, model->centroid);
		vecmat(model->latmat,r);
		ARR3SUB(r, point);
		r2 = VEC3MAGSQ(r);
		
		/*  ARR3SET(coreCartesian, core->x);
		vecmat(model->latmat,coreCartesian);
		ARR3SUB(coreCartesian,x);
		r2 = VEC3MAGSQ(r);
		*/  
		rhoNugget = hfs_calc_density(core->atom_code, r2);
		rhoMol += rhoNugget;
	}


#ifdef JJM_DEBUG02
	if ( rhoMol > 0.0 ){
		fprintf(stderr,"rhoMol: %.3f\n",rhoMol);
	}
#endif
	
	
	return(rhoMol);
}
