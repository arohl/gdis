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
/**************/
/* structures */
/**************/

struct surface_property_pak
{
	float	de;				/* distance to the nearest exernal atom */
	int		eType;			/* atomic number of the nearest external atom */
	float	di;				/* distance to the nearest internal atom */
	int		iType;			/* atomic number of the nearest internal atom */
	float	shape_index;	
	float	curvedness;		/* properties based on surface curvature */
	float	promol_den;	/* promolecule density, for Hirshfeld surface (== half pro crystal density) */
};

struct smv_pak
{
    gdouble value;
    gdouble colour[3];
    
    /* Added by JJM */
    gfloat property;						/* property should eventually be replaced by the struct */
	struct surface_property_pak properties;
    
    /* coords */
    gdouble x[3];
    gdouble rx[3];
    
    /* normal */
    gdouble n[3];
    gdouble nx[3];
    
    /* adjacent point */
    GSList *adj;
};

struct smt_pak
{
    struct smv_pak *point[3];
};



/**************/
/* prototypes */
/**************/

//void ms_properties(gint method, GSList *ms_points, struct model_pak *model);
void ms_colour_surface(gint, gint, GSList *, struct model_pak *);

void ms_dock_colour(gdouble *, gdouble, gdouble, gdouble);
void ms_epot_colour(gdouble *, gdouble, gdouble, gdouble);
void ms_afm_colour(gdouble *, gdouble, struct model_pak *);
void ms_hfs_colour(gdouble *, gdouble);
void ms_colour_by_de(GSList *, struct model_pak *);

void ms_colour_by_shape(GSList *points, struct model_pak *model, gint propertyType, gint method);


void ms_cube(gdouble, gint, gint, struct model_pak *);



/****************************************************************/
/* I've included these here because I am using them in molsurf, */
/* but they should probably go somewhere more sensible (JJM)	*/
/****************************************************************/
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/****************************************************************/
/* Enumerated types for surface properties						*/
/*																*/ 
/* These values currently determine which property is 			*/
/* calculated, but they will eventually flag which property		*/
/* is currently displayed.										*/
/****************************************************************/
enum {	CURVEDNESS,			/* */
		SHAPE_INDEX,
		DE,					/* distance to the nearest atom outside the surface */
		DI,					/* distance to the nearest atom inside the surface */
		ETYPE,				/* atomic number of the nearest atom outside the surface */
		ITYPE,				/* atomic number of the nearest atom inside the surface */
		PROMOL_DEN,			/* total electron density from the promolecule (equals half the total procrystal density */ 
		AFM,				
		ELPOT,
		};
