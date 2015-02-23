/*
 *  Brute force symmetry analyzer.
 *  This is actually C++ program, masquerading as a C one!
 *
 *  (C) 1996, 2003 S. Patchkovskii, Serguei.Patchkovskii@sympatico.ca
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Modified by Sean Fleming to use glib, and for inclusion in gdis */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gdis.h"
#include "coords.h"
#include "edit.h"
#include "matrix.h"
#include "sginfo.h"
#include "gui_shorts.h"
#include "interface.h"
#include "dialog.h"
#include "opengl.h"

#define	DIMENSION 3
#define MAXPARAM  7
#define DEBUG 0

extern struct sysenv_pak sysenv;

GtkListStore *symmetry_ls=NULL;
GtkWidget *symmetry_tv=NULL;

typedef struct 
  {
  gint     type;
  gdouble  x[DIMENSION];
  } OBJECT;

/*
 *  All specific structures should have corresponding elements in the
 *  same position generic structure does.
 *
 *  Planes are characterized by the surface normal direction
 *  (taken in the direction *from* the coordinate origin)
 *  and distance from the coordinate origin to the plane
 *  in the direction of the surface normal.
 *
 *  Inversion is characterized by location of the inversion center.
 *
 *  Rotation is characterized by a vector (distance+direction) from the origin
 *  to the rotation axis, axis direction and rotation order. Rotations
 *  are in the clockwise direction looking opposite to the direction
 *  of the axis. Note that this definition of the rotation axis
 *  is *not* unique, since an arbitrary multiple of the axis direction
 *  can be added to the position vector without changing actual operation.
 *
 *  Mirror rotation is defined by the same parameters as normal rotation,
 *  but the origin is now unambiguous since it defines the position of the
 *  plane associated with the axis.
 *
 */

typedef struct _SYMMETRY_ELEMENT_ 
  {
  void (*transform_atom) 
       (struct _SYMMETRY_ELEMENT_ *el, OBJECT *from, OBJECT *to );
  gint *transform;   /* Correspondence table for the transformation  */
  gint order;        /* Applying transformation this many times is identity */
  gint nparam;       /* 4 for inversion and planes, 7 for axes */
  gdouble maxdev;    /* Larges error associated with the element */
  gdouble distance ;
  gdouble normal[DIMENSION];
  gdouble direction[DIMENSION];
  } SYMMETRY_ELEMENT;

typedef struct 
  {
  gchar *group_name;        /* Canonical group name */
  gchar *symmetry_code;     /* Group symmetry code  */
  gint (*check)(void);      /* Additional verification routine, not used */
  } POINT_GROUP;

/****************/
/* globals, yuk */
/****************/
gdouble                 ToleranceSame         = 1e-3 ;
gdouble                 TolerancePrimary      = 5e-2 ;
gdouble                 ToleranceFinal        = 1e-4 ;
gdouble                 MaxOptStep            = 5e-1 ;
gdouble                 MinOptStep            = 1e-7 ;
gdouble                 GradientStep          = 1e-7 ;
gdouble                 OptChangeThreshold    = 1e-10 ;
gdouble                 CenterOfSomething[ DIMENSION ] ;
gdouble *               DistanceFromCenter    = NULL ;
gint                    verbose               = 0 ;
gint                    MaxOptCycles          = 200 ;
gint                    OptChangeHits         = 5 ;
gint                    MaxAxisOrder          = 20 ;
gint                    AtomsCount            = 0 ;
OBJECT *                 Atoms                 = NULL ;
gint                    PlanesCount           = 0 ;
SYMMETRY_ELEMENT **    Planes                = NULL ;
SYMMETRY_ELEMENT *     MolecularPlane        = NULL ;
gint                    InversionCentersCount = 0 ;
SYMMETRY_ELEMENT **    InversionCenters      = NULL ;
gint                    NormalAxesCount       = 0 ;
SYMMETRY_ELEMENT **    NormalAxes            = NULL ;
gint                    ImproperAxesCount     = 0 ;
SYMMETRY_ELEMENT **    ImproperAxes          = NULL ;
gint *                  NormalAxesCounts      = NULL ;
gint *                  ImproperAxesCounts    = NULL ;
gint                    BadOptimization       = 0 ;
gchar *                 SymmetryCode          = NULL ;
/*
 *    Statistics
 */
gint                   StatTotal             = 0 ;
gint                   StatEarly             = 0 ;
gint                   StatPairs             = 0 ;
gint                   StatDups              = 0 ;
gint                   StatOrder             = 0 ;
gint                   StatOpt               = 0 ;
gint                   StatAccept            = 0 ;

/*
 *    Point groups I know about
 */
gint true(void){return 1;}

POINT_GROUP PointGroups[] = 
  {
  {"C1",    "",                                                       true},
  {"Cs",    "(sigma) ",                                               true},
  {"Ci",    "(i) ",                                                   true},
  {"C2",    "(C2) ",                                                  true},
  {"C3",    "(C3) ",                                                  true},
  {"C4",    "(C4) (C2) ",                                             true},
  {"C5",    "(C5) ",                                                  true},
  {"C6",    "(C6) (C3) (C2) ",                                        true},
  {"C7",    "(C7) ",                                                  true},
  {"C8",    "(C8) (C4) (C2) ",                                        true},
  {"D2",    "3*(C2) ",                                                true},
  {"D3",    "(C3) 3*(C2) ",                                           true},
  {"D4",    "(C4) 5*(C2) ",                                           true},
  {"D5",    "(C5) 5*(C2) ",                                           true},
  {"D6",    "(C6) (C3) 7*(C2) ",                                      true},
  {"D7",    "(C7) 7*(C2) ",                                           true},
  {"D8",    "(C8) (C4) 9*(C2) ",                                      true},
  {"C2v",   "(C2) 2*(sigma) ",                                        true},
  {"C3v",   "(C3) 3*(sigma) ",                                        true},
  {"C4v",   "(C4) (C2) 4*(sigma) ",                                   true},
  {"C5v",   "(C5) 5*(sigma) ",                                        true},
  {"C6v",   "(C6) (C3) (C2) 6*(sigma) ",                              true},
  {"C7v",   "(C7) 7*(sigma) ",                                        true},
  {"C8v",   "(C8) (C4) (C2) 8*(sigma) ",                              true},
  {"C2h",   "(i) (C2) (sigma) ",                                      true},
  {"C3h",   "(C3) (S3) (sigma) ",                                     true},
  {"C4h",   "(i) (C4) (C2) (S4) (sigma) ",                            true},
  {"C5h",   "(C5) (S5) (sigma) ",                                     true},
  {"C6h",   "(i) (C6) (C3) (C2) (S6) (S3) (sigma) ",                  true},
  {"C7h",   "(C7) (S7) (sigma) ",                                     true},
  {"C8h",   "(i) (C8) (C4) (C2) (S8) (S4) (sigma) ",                  true},
  {"D2h",   "(i) 3*(C2) 3*(sigma) ",                                  true},
  {"D3h",   "(C3) 3*(C2) (S3) 4*(sigma) ",                            true},
  {"D4h",   "(i) (C4) 5*(C2) (S4) 5*(sigma) ",                        true},
  {"D5h",   "(C5) 5*(C2) (S5) 6*(sigma) ",                            true},
  {"D6h",   "(i) (C6) (C3) 7*(C2) (S6) (S3) 7*(sigma) ",              true},
  {"D7h",   "(C7) 7*(C2) (S7) 8*(sigma) ",                            true},
  {"D8h",   "(i) (C8) (C4) 9*(C2) (S8) (S4) 9*(sigma) ",              true},
  {"D2d",   "3*(C2) (S4) 2*(sigma) ",                                 true},
  {"D3d",   "(i) (C3) 3*(C2) (S6) 3*(sigma) ",                        true},
  {"D4d",   "(C4) 5*(C2) (S8) 4*(sigma) ",                            true},
  {"D5d",   "(i) (C5) 5*(C2) (S10) 5*(sigma) ",                       true},
  {"D6d",   "(C6) (C3) 7*(C2) (S12) (S4) 6*(sigma) ",                 true},
  {"D7d",   "(i) (C7) 7*(C2) (S14) 7*(sigma) ",                       true},
  {"D8d",   "(C8) (C4) 9*(C2) (S16) 8*(sigma) ",                      true},
  {"S4",    "(C2) (S4) ",                                             true},
  {"S6",    "(i) (C3) (S6) ",                                         true},
  {"S8",    "(C4) (C2) (S8) ",                                        true},
  {"T",     "4*(C3) 3*(C2) ",                                         true},
  {"Th",    "(i) 4*(C3) 3*(C2) 4*(S6) 3*(sigma) ",                    true},
  {"Td",    "4*(C3) 3*(C2) 3*(S4) 6*(sigma) ",                        true},
  {"O",     "3*(C4) 4*(C3) 9*(C2) ",                                  true},
  {"Oh",    "(i) 3*(C4) 4*(C3) 9*(C2) 4*(S6) 3*(S4) 9*(sigma) ",      true},
  {"Cinfv", "(Cinf) (sigma) ",                                        true},
  {"Dinfh", "(i) (Cinf) (C2) 2*(sigma) ",                             true},
  {"I",     "6*(C5) 10*(C3) 15*(C2) ",                                true},
  {"Ih",    "(i) 6*(C5) 10*(C3) 15*(C2) 6*(S10) 10*(S6) 15*(sigma) ", true},
  {"Kh",    "(i) (Cinf) (sigma) ",                                    true},
  };


#define PointGroupsCount (sizeof(PointGroups)/sizeof(POINT_GROUP))
gchar *                 PointGroupRejectionReason = NULL ;

/*
 *   Generic functions
 */

gdouble pow2( gdouble x )
{
return x * x ;
}

gint establish_pairs( SYMMETRY_ELEMENT *elem )
{
gint i, j, k, best_j ;
gchar *atom_used = g_malloc0(AtomsCount);
gdouble distance, best_distance;
OBJECT symmetric;

for( i = 0 ; i < AtomsCount ; i++ ){
    if( elem->transform[i] >= AtomsCount ){ /* No symmetric atom yet          */
        if( verbose > 2 ) printf( "        looking for a pair for %d\n", i ) ;
        elem->transform_atom( elem, Atoms+i, &symmetric ) ;
        if( verbose > 2 ) printf( "        new coordinates are: (%g,%g,%g)\n", 
                              symmetric.x[0], symmetric.x[1], symmetric.x[2] ) ;
        best_j        = i ;
        best_distance = 2*TolerancePrimary ;/* Performance value we'll reject */
        for( j = 0 ; j < AtomsCount ; j++ ){
            if( Atoms[j].type != symmetric.type || atom_used[j] )
                continue ;
            for( k = 0, distance = 0 ; k < DIMENSION ; k++ ){
                distance += pow2( symmetric.x[k] - Atoms[j].x[k] ) ;
                }
            distance = sqrt( distance ) ;
            if( verbose > 2 ) printf( "        distance to %d is %g\n", j, distance ) ;
            if( distance < best_distance ){
                best_j        = j ;
                best_distance = distance ;
                }
            }
        if( best_distance > TolerancePrimary ){ /* Too bad, there is no symmetric atom */
            if( verbose > 0 ) 
                printf( "        no pair for atom %d - best was %d with err = %g\n", i, best_j, best_distance ) ;
            g_free( atom_used ) ;
            return -1 ;
            }
        elem->transform[i] = best_j ;
        atom_used[best_j]  = 1 ;
        if( verbose > 1 ) printf( "        atom %d transforms to the atom %d, err = %g\n", i, best_j, best_distance ) ;
        }
    }
g_free( atom_used ) ;
return 0 ;
}

gint check_transform_order( SYMMETRY_ELEMENT *elem )
{
gint             i, j, k ;
void            rotate_reflect_atom( SYMMETRY_ELEMENT *, OBJECT *, OBJECT *) ;

for( i = 0 ; i < AtomsCount ; i++ ){
    if( elem->transform[i] == i )   /* Identity transform is Ok for any order */
        continue ;
    if( elem->transform_atom == rotate_reflect_atom ){
        j = elem->transform[i] ;
        if( elem->transform[j] == i )
            continue ; /* Second-order transform is Ok for improper axis */
        }
    for( j = elem->order - 1, k = elem->transform[i] ; j > 0 ; j--, k = elem->transform[k] ){
        if( k == i ){
            if( verbose > 0 ) printf( "        transform looped %d steps too early from atom %d\n", j, i ) ;
            return -1 ;
            }
        }
    if( k != i && elem->transform_atom == rotate_reflect_atom ){
        /* For improper axes, the complete loop may also take twice the order */
        for( j = elem->order ; j > 0 ; j--, k = elem->transform[k] ){
            if( k == i ){
                if( verbose > 0 ) printf( "        (improper) transform looped %d steps too early from atom %d\n", j, i ) ;
                return -1 ;
                }
            }
        }
    if( k != i ){
        if( verbose > 0 ) printf( "        transform failed to loop after %d steps from atom %d\n", elem->order, i ) ;
        return -1 ;
        }
    }
return 0 ;
}

gint same_transform( SYMMETRY_ELEMENT *a, SYMMETRY_ELEMENT *b )
{
gint               i, j ;
gint               code ;

if( ( a->order != b->order ) || ( a->nparam != b->nparam ) || ( a->transform_atom != b->transform_atom ) )
    return 0 ;
for( i = 0, code = 1 ; i < AtomsCount ; i++ ){
    if( a->transform[i] != b->transform[i] ){
        code = 0 ;
        break ;
        }
    }
if( code == 0 && a->order > 2 ){  /* b can also be a reverse transformation for a */
    for( i = 0 ; i < AtomsCount ; i++ ){
        j = a->transform[i] ;
        if( b->transform[j] != i )
            return 0 ;
        }
    return 1 ;
    }
return code ;
}

SYMMETRY_ELEMENT *alloc_symmetry_element( void )
{
gint                i;
SYMMETRY_ELEMENT *elem;

elem = g_malloc0(sizeof(SYMMETRY_ELEMENT));
elem->transform = g_malloc0(AtomsCount*sizeof(int));

for (i=0 ; i<AtomsCount ; i++)
  elem->transform[i] = AtomsCount + 1 ; /* An impossible value */

return elem ;
}

void destroy_symmetry_element( SYMMETRY_ELEMENT *elem )
{
if (elem != NULL)
  {
  if (elem->transform != NULL)
    g_free(elem->transform);
  g_free(elem);
  }
}

gint check_transform_quality( SYMMETRY_ELEMENT *elem )
{
gint               i, j, k ;
OBJECT              symmetric ;
gdouble            r, max_r ;

for( i = 0, max_r = 0 ; i < AtomsCount ; i++ ){
    j = elem->transform[i] ;
    elem->transform_atom( elem, Atoms + i, &symmetric ) ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( symmetric.x[k] - Atoms[j].x[k] ) ;
        }
    r = sqrt( r ) ;
    if( r > ToleranceFinal ){
        if( verbose > 0 ) printf( "        distance to symmetric atom (%g) is too big for %d\n", r, i ) ;
        return -1 ;
        }
    if( r > max_r ) max_r = r ;
    }
elem->maxdev = max_r ;
return 0 ;
}

gdouble eval_optimization_target_function(SYMMETRY_ELEMENT *elem, gint *finish)
{
gint               i, j, k ;
OBJECT              symmetric ;
gdouble            target, r, maxr ;

if( elem->nparam >= 4 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->normal[k]*elem->normal[k] ;
        }
    r = sqrt( r ) ;
    if( r < ToleranceSame ){
        fprintf( stderr, "Normal collapsed!\n" ) ;
        return(0.0);
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->normal[k] /= r ;
        }
    if( elem->distance < 0 ){
        elem->distance = -elem->distance ;
        for( k = 0 ; k < DIMENSION ; k++ ){
            elem->normal[k] = -elem->normal[k] ;
            }
        }
    }
if( elem->nparam >= 7 ){
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += elem->direction[k]*elem->direction[k] ;
        }
    r = sqrt( r ) ;
    if( r < ToleranceSame ){
        fprintf( stderr, "Direction collapsed!\n" ) ;
        return(0.0);
        }
    for( k = 0 ; k < DIMENSION ; k++ ){
        elem->direction[k] /= r ;
        }
    }
for( i = 0, target = maxr = 0 ; i < AtomsCount ; i++ ){
    elem->transform_atom( elem, Atoms + i, &symmetric ) ;
    j = elem->transform[i] ;
    for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
        r += pow2( Atoms[j].x[k] - symmetric.x[k] ) ;
        }
    if( r > maxr ) maxr = r ;
    target += r ;
    }
if( finish != NULL ){
    *finish = 0 ;
    if( sqrt( maxr ) < ToleranceFinal )
        *finish = 1 ;
    }
return target ;
}

void get_params( SYMMETRY_ELEMENT *elem, gdouble values[] )
{
memcpy( values, &elem->distance, elem->nparam * sizeof( gdouble ) ) ;
}

void set_params( SYMMETRY_ELEMENT *elem, gdouble values[] )
{
memcpy( &elem->distance, values, elem->nparam * sizeof( gdouble ) ) ;
}

gint optimize_transformation_params(SYMMETRY_ELEMENT *elem)
{
gdouble            values[ MAXPARAM ] ;
gdouble            grad  [ MAXPARAM ] ;
gdouble            force [ MAXPARAM ] ;
gdouble            step  [ MAXPARAM ] ;
gdouble            f, fold, fnew, fnew2, fdn, fup, snorm ;
gdouble            a, b, x ;
gint               vars  = elem->nparam ;
gint               cycle = 0 ;
gint               i, finish ;
gint               hits = 0 ;

if (vars > MAXPARAM)
  {
  fprintf(stderr, "Catastrophe in optimize_transformation_params()!\n");
  return(1);
  }

f = 0 ;
do {
    fold = f ;
    f    = eval_optimization_target_function( elem, &finish ) ;
    /* Evaluate function, gradient and diagonal force constants */
    if( verbose > 1 ) printf( "            function value = %g\n", f ) ;
    if( finish ){
        if( verbose > 1 ) printf( "        function value is small enough\n" ) ;
        break ;
        }
    if( cycle > 0 ){
        if( fabs( f-fold ) > OptChangeThreshold )
             hits = 0 ;
        else hits++ ;
        if( hits >= OptChangeHits ){
            if( verbose > 1 ) printf( "        no progress is made, stop optimization\n" ) ;
            break ;
            }
        }
    get_params( elem, values ) ;
    for( i = 0 ; i < vars ; i++ ){
        values[i] -= GradientStep ;
        set_params( elem, values ) ;
        fdn        = eval_optimization_target_function( elem, NULL ) ;
        values[i] += 2*GradientStep ;
        set_params( elem, values ) ;
        fup        = eval_optimization_target_function( elem, NULL ) ;
        values[i] -= GradientStep ;
        grad[i]    = ( fup - fdn ) / ( 2 * GradientStep ) ;
        force[i]   = ( fup + fdn - 2*f ) / ( GradientStep * GradientStep ) ;
        if( verbose > 1 ) printf( "        i = %d, grad = %12.6e, force = %12.6e\n", i, grad[i], force[i] ) ;
        }
    /* Do a quasy-Newton step */
    for( i = 0, snorm = 0 ; i < vars ; i++ ){
        if( force[i] <  0   ) force[i] = -force[i] ;
        if( force[i] < 1e-3 ) force[i] = 1e-3 ;
        if( force[i] > 1e3  ) force[i] = 1e3 ;
        step[i] = - grad[i]/force[i] ;
        snorm += step[i] * step[i] ;
        }
    snorm = sqrt( snorm ) ;
    if( snorm > MaxOptStep ){ /* Renormalize step */
        for( i = 0 ; i < vars ; i++ )
            step[i] *= MaxOptStep/snorm ;
        snorm = MaxOptStep ;
        }
    do {
        for( i = 0 ; i < vars ; i++ ){
            values[i] += step[i] ;
            }
        set_params( elem, values ) ;
        fnew = eval_optimization_target_function( elem, NULL ) ;
        if( fnew < f )
            break ;
        for( i = 0 ; i < vars ; i++ ){
            values[i] -= step[i] ;
            step  [i] /= 2 ;
            }
        set_params( elem, values ) ;
        snorm /= 2 ;
        } while( snorm > MinOptStep ) ;
        if( (snorm > MinOptStep) && (snorm < MaxOptStep / 2) ){  /* try to do quadratic interpolation */
            for( i = 0 ; i < vars ; i++ )
                values[i] += step[i] ;
            set_params( elem, values ) ;
            fnew2 = eval_optimization_target_function( elem, NULL ) ;
            if( verbose > 1 ) printf( "        interpolation base points: %g, %g, %g\n", f, fnew, fnew2 ) ;
            for( i = 0 ; i < vars ; i++ )
                values[i] -= 2*step[i] ;
            a     = ( 4*f - fnew2 - 3*fnew ) / 2 ;
            b     = ( f + fnew2 - 2*fnew ) / 2 ;
            if( verbose > 1 ) printf( "        linear interpolation coefficients %g, %g\n", a, b ) ;
            if( b > 0 ){
                x = -a/(2*b) ;
                if( x > 0.2 && x < 1.8 ){
                    if( verbose > 1 ) printf( "        interpolated: %g\n", x ) ;
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += x*step[i] ;
                    }
                else b = 0 ;
                }
            if( b <= 0 ){
                if( fnew2 < fnew ){
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += 2*step[i] ;
                    }
                else {
                    for( i = 0 ; i < vars ; i++ )
                        values[i] += step[i] ;
                    }
                }
            set_params( elem, values ) ;
            }
    } while( snorm > MinOptStep && ++cycle < MaxOptCycles ) ;
f = eval_optimization_target_function( elem, NULL ) ;
if( cycle >= MaxOptCycles ) BadOptimization = 1 ;
if( verbose > 0 ) {
    if( cycle >= MaxOptCycles )
        printf( "        maximum number of optimization cycles made\n" ) ;
        printf( "        optimization completed after %d cycles with f = %g\n", cycle, f ) ;
    }
return(0);
}

gint refine_symmetry_element( SYMMETRY_ELEMENT *elem, gint build_table )
{
        gint               i ;


if( build_table && (establish_pairs( elem ) < 0) ){
    StatPairs++ ;
    if( verbose > 0 ) printf( "        no transformation correspondence table can be constructed\n" ) ;
    return -1 ;
    }
for( i = 0 ; i < PlanesCount ; i++ ){
    if( same_transform( Planes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to plane %d\n", i ) ;
        return -1 ;
        }
    }

for( i = 0 ; i < InversionCentersCount ; i++ ){
    if( same_transform( InversionCenters[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to inversion center %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < NormalAxesCount ; i++ ){
    if( same_transform( NormalAxes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to normal axis %d\n", i ) ;
        return -1 ;
        }
    }
for( i = 0 ; i < ImproperAxesCount ; i++ ){
    if( same_transform( ImproperAxes[i], elem ) ){
        StatDups++ ;
        if( verbose > 0 ) printf( "        transformation is identical to improper axis %d\n", i ) ;
        return -1 ;
        }
    }
if( check_transform_order( elem ) < 0 ){
    StatOrder++ ;
    if( verbose > 0 ) printf( "        incorrect transformation order\n" ) ;
    return -1 ;
    }

if (optimize_transformation_params( elem ))
  return(-1);

if( check_transform_quality( elem ) < 0 ){
    StatOpt++ ;
    if( verbose > 0 ) printf( "        refined transformation does not pass the numeric threshold\n" ) ;
    return -1 ;
    }
StatAccept++ ;
return 0 ;
}

/*
 *   Plane-specific functions
 */

void mirror_atom( SYMMETRY_ELEMENT *plane, OBJECT *from, OBJECT *to )
{
gint                i ;
gdouble             r ;

for( i = 0, r = plane->distance ; i < DIMENSION ; i++ ){
    r -= from->x[i] * plane->normal[i] ;
    }
to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = from->x[i] + 2*r*plane->normal[i] ;
    }
}

SYMMETRY_ELEMENT * init_mirror_plane( gint i, gint j )
{
SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
gdouble             dx[ DIMENSION ], midpoint[ DIMENSION ], rab, r ;
gint                k ;

if( verbose > 0 ) printf( "Trying mirror plane for atoms %d,%d\n", i, j ) ;
StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 2 ;
plane->nparam         = 4 ;
for( k = 0, rab = 0 ; k < DIMENSION ; k++ ){
    dx[k]       = Atoms[i].x[k] - Atoms[j].x[k] ;
    midpoint[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2.0 ;
    rab        += dx[k]*dx[k] ;
    }
rab = sqrt(rab) ;
if (rab < ToleranceSame)
  {
  fprintf(stderr, "Atoms %d and %d coincide (r = %g)\n", i, j, rab );
  return NULL;
  }

for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    plane->normal[k] = dx[k]/rab ;
    r += midpoint[k]*plane->normal[k] ;
    }
if( r < 0 ){  /* Reverce normal direction, distance is always positive! */
    r = -r ;
    for( k = 0 ; k < DIMENSION ; k++ ){
        plane->normal[k] = -plane->normal[k] ;
        }
    }
plane->distance = r ;
if( verbose > 0 ) printf( "    initial plane is at %g from the origin\n", r ) ;
if( refine_symmetry_element( plane, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}

SYMMETRY_ELEMENT *init_ultimate_plane( void )
{
SYMMETRY_ELEMENT * plane = alloc_symmetry_element() ;
gdouble             d0[ DIMENSION ], d1[ DIMENSION ], d2[ DIMENSION ] ;
gdouble             p[ DIMENSION ] ;
gdouble             r, s0, s1, s2 ;
gdouble *           d ;
gint                i, j, k ;

if( verbose > 0 ) printf( "Trying whole-molecule mirror plane\n" ) ;
StatTotal++ ;
plane->transform_atom = mirror_atom ;
plane->order          = 1 ;
plane->nparam         = 4 ;
for( k = 0 ; k < DIMENSION ; k++ )
    d0[k] = d1[k] = d2[k] = 0 ;
d0[0] = 1 ; d1[1] = 1 ; d2[2] = 1 ;
for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            p[k] = Atoms[i].x[k] - Atoms[j].x[k] ;
            r   += p[k]*p[k] ;
            }
        r = sqrt(r) ;
        for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
            p[k] /= r ;
            s0   += p[k]*d0[k] ;
            s1   += p[k]*d1[k] ;
            s2   += p[k]*d2[k] ;
            }
        for( k = 0 ; k < DIMENSION ; k++ ){
            d0[k] -= s0*p[k] ;
            d1[k] -= s1*p[k] ;
            d2[k] -= s2*p[k] ;
            }
        }
    }
for( k = 0, s0=s1=s2=0 ; k < DIMENSION ; k++ ){
    s0 += d0[k] ;
    s1 += d1[k] ;
    s2 += d2[k] ;
    }
d = NULL ;
if( s0 >= s1 && s0 >= s2 ) d = d0 ;
if( s1 >= s0 && s1 >= s2 ) d = d1 ;
if( s2 >= s0 && s2 >= s1 ) d = d2 ;

if(d == NULL)
  {
  fprintf(stderr, "ERROR: init_ultimate_plane(): %g, %g and %g have no ordering!\n", s0, s1, s2);
  return NULL;
  }

for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += d[k]*d[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        plane->normal[k] = d[k]/r ;
    }
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        plane->normal[k] = 0 ;
    plane->normal[0] = 1 ;
    }
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*plane->normal[k] ;
plane->distance = r ;
for( k = 0 ; k < AtomsCount ; k++ )
    plane->transform[k] = k ;
if( refine_symmetry_element( plane, 0 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the plane\n" ) ;
    destroy_symmetry_element( plane ) ;
    return NULL ;
    }
return plane ;
}

/*
 *   Inversion-center specific functions
 */
void invert_atom( SYMMETRY_ELEMENT *center, OBJECT *from, OBJECT *to )
{
gint                i ;

to->type = from->type ;
for( i = 0 ; i < DIMENSION ; i++ ){
    to->x[i] = 2*center->distance*center->normal[i] - from->x[i] ;
    }
}

SYMMETRY_ELEMENT *init_inversion_center( void )
{
SYMMETRY_ELEMENT * center = alloc_symmetry_element() ;
gint                k ;
gdouble             r ;

if( verbose > 0 ) printf( "Trying inversion center at the center of something\n" ) ;
StatTotal++ ;
center->transform_atom = invert_atom ;
center->order          = 2 ;
center->nparam         = 4 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        center->normal[k] = CenterOfSomething[k]/r ;
    }
else {
    center->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        center->normal[k] = 0 ;
    }
center->distance = r ;
if( verbose > 0 ) printf( "    initial inversion center is at %g from the origin\n", r ) ;
if( refine_symmetry_element( center, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the inversion center\n" ) ;
    destroy_symmetry_element( center ) ;
    return NULL ;
    }
return center ;
}

/*
 *   Normal rotation axis-specific routines.
 */
void rotate_atom(SYMMETRY_ELEMENT *axis, OBJECT *from, OBJECT *to)
{
gdouble             x[3], y[3], a[3], b[3], c[3] ;
gdouble             angle = axis->order ? 2*PI/axis->order : 1.0 ;
gdouble             a_sin = sin( angle ) ;
gdouble             a_cos = cos( angle ) ;
gdouble             dot ;
gint                i ;

if (DIMENSION != 3)
  {
  fprintf(stderr, "Catastrophe in rotate_atom!\n");
  return;
  }

for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type;
}

SYMMETRY_ELEMENT *init_ultimate_axis(void)
{
SYMMETRY_ELEMENT * axis = alloc_symmetry_element() ;
gdouble             dir[ DIMENSION ], rel[ DIMENSION ] ;
gdouble             s ;
gint                i, k ;

if( verbose > 0 ) printf( "Trying infinity axis\n" ) ;
StatTotal++ ;
axis->transform_atom = rotate_atom ;
axis->order          = 0 ;
axis->nparam         = 7 ;
for( k = 0 ; k < DIMENSION ; k++ )
    dir[k] = 0 ;
for( i = 0 ; i < AtomsCount ; i++ ){
    for( k = 0, s = 0 ; k < DIMENSION ; k++ ){
        rel[k] = Atoms[i].x[k] - CenterOfSomething[k] ;
        s     += rel[k]*dir[k] ;
        }
    if( s >= 0 )
         for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] += rel[k] ;
    else for( k = 0 ; k < DIMENSION ; k++ )
             dir[k] -= rel[k] ;
    }
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( dir[k] ) ;
s = sqrt(s) ;
if( s > 0 )
     for( k = 0 ; k < DIMENSION ; k++ )
         dir[k] /= s ;
else dir[0] = 1 ;
for( k = 0 ; k < DIMENSION ; k++ )
    axis->direction[k] = dir[k] ;
for( k = 0, s = 0 ; k < DIMENSION ; k++ )
    s += pow2( CenterOfSomething[k] ) ;
s = sqrt(s) ;
if( s > 0 )
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = CenterOfSomething[k]/s ;
else {
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    axis->normal[0] = 1 ;
    }
axis->distance = s ;
for( k = 0 ; k < AtomsCount ; k++ )
    axis->transform[k] = k ;
if( refine_symmetry_element( axis, 0 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the infinity axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}


SYMMETRY_ELEMENT *init_c2_axis( gint i, gint j, gdouble support[ DIMENSION ] )
{
SYMMETRY_ELEMENT * axis ;
gint                k ;
gdouble             ris, rjs ;
gdouble             r, center[ DIMENSION ] ;

if( verbose > 0 ) 
    printf( "Trying c2 axis for the pair (%d,%d) with the support (%g,%g,%g)\n", 
             i, j, support[0], support[1], support[2] ) ;
StatTotal++ ;
/* First, do a quick sanity check */
for( k = 0, ris = rjs = 0 ; k < DIMENSION ; k++ ){
    ris += pow2( Atoms[i].x[k] - support[k] ) ;
    rjs += pow2( Atoms[j].x[k] - support[k] ) ;
    }
ris = sqrt( ris ) ;
rjs = sqrt( rjs ) ;
if( fabs( ris - rjs ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    Support can't actually define a rotation axis\n" ) ;
    return NULL ;
    }
axis                 = alloc_symmetry_element() ;
axis->transform_atom = rotate_atom ;
axis->order          = 2 ;
axis->nparam         = 7 ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ )
    r += CenterOfSomething[k]*CenterOfSomething[k] ;
r = sqrt(r) ;
if( r > 0 ){
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->normal[k] = CenterOfSomething[k]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( k = 1 ; k < DIMENSION ; k++ )
        axis->normal[k] = 0 ;
    }
axis->distance = r ;
for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
    center[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2 - support[k] ;
    r        += center[k]*center[k] ;
    }
r = sqrt(r) ;
if( r <= TolerancePrimary ){ /* c2 is underdefined, let's do something special */
    if( MolecularPlane != NULL ){
        if( verbose > 0 ) printf( "    c2 is underdefined, but there is a molecular plane\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] = MolecularPlane->normal[k] ;
        }
    else {
        if( verbose > 0 ) printf( "    c2 is underdefined, trying random direction\n" ) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            center[k] = Atoms[i].x[k] - Atoms[j].x[k] ;
        if( fabs( center[2] ) + fabs( center[1] ) > ToleranceSame ){
            axis->direction[0] =  0 ;
            axis->direction[1] =  center[2] ;
            axis->direction[2] = -center[1] ;
            }
        else {
            axis->direction[0] = -center[2] ;
            axis->direction[1] =  0 ;
            axis->direction[2] =  center[0] ;
            }
        for( k = 0, r = 0 ; k < DIMENSION ; k++ )
            r += axis->direction[k] * axis->direction[k] ;
        r = sqrt(r) ;
        for( k = 0 ; k < DIMENSION ; k++ )
            axis->direction[k] /= r ;
        }
    }
else { /* direction is Ok, renormalize it */
    for( k = 0 ; k < DIMENSION ; k++ )
        axis->direction[k] = center[k]/r ;
    }
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the c2 axis\n" ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

SYMMETRY_ELEMENT *init_axis_parameters
                  (gdouble a[3], gdouble b[3], gdouble c[3])
{
        SYMMETRY_ELEMENT * axis ;
        gint                i, order, sign ;
        gdouble             ra, rb, rc, rab, rbc, rac, r ;
        gdouble             angle ;

ra = rb = rc = rab = rbc = rac = 0 ;
for( i = 0 ; i < DIMENSION ; i++ ){
    ra  += a[i]*a[i] ;
    rb  += b[i]*b[i] ;
    rc  += c[i]*c[i] ;
    }
ra = sqrt(ra) ; rb  = sqrt(rb) ; rc  = sqrt(rc) ;
if( fabs( ra - rb ) > TolerancePrimary || fabs( ra - rc ) > TolerancePrimary || fabs( rb - rc ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    points are not on a sphere\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ ){
    rab += (a[i]-b[i])*(a[i]-b[i]) ;
    rac += (a[i]-c[i])*(a[i]-c[i]) ;
    rbc += (c[i]-b[i])*(c[i]-b[i]) ;
    }
rab = sqrt(rab) ;
rac = sqrt(rac) ;
rbc = sqrt(rbc) ;
if( fabs( rab - rbc ) > TolerancePrimary ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    points can't be rotation-equivalent\n" ) ;
    return NULL ;
    }
if( rab <= ToleranceSame || rbc <= ToleranceSame || rac <= ToleranceSame ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    rotation is underdefined by these points\n" ) ;
    return NULL ;
    }
rab   = (rab+rbc)/2 ;
angle = PI - 2*asin( rac/(2*rab) ) ;
if( verbose > 1 ) printf( "    rotation angle is %f\n", angle ) ;
if( fabs(angle) <= PI/(MaxAxisOrder+1) ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    atoms are too close to a straight line\n" ) ;
    return NULL ;
    }
order = floor( (2*PI)/angle + 0.5 ) ;
if( order <= 2 || order > MaxAxisOrder ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    rotation axis order (%d) is not from 3 to %d\n", order, MaxAxisOrder ) ;
    return NULL ;
    }
axis = alloc_symmetry_element() ;
axis->order          = order ;
axis->nparam         = 7 ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += CenterOfSomething[i]*CenterOfSomething[i] ;
r = sqrt(r) ;
if( r > 0 ){
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->normal[i] = CenterOfSomething[i]/r ;
    }
else {
    axis->normal[0] = 1 ;
    for( i = 1 ; i < DIMENSION ; i++ )
        axis->normal[i] = 0 ;
    }
axis->distance = r ;
axis->direction[0] = (b[1]-a[1])*(c[2]-b[2]) - (b[2]-a[2])*(c[1]-b[1]) ;
axis->direction[1] = (b[2]-a[2])*(c[0]-b[0]) - (b[0]-a[0])*(c[2]-b[2]) ;
axis->direction[2] = (b[0]-a[0])*(c[1]-b[1]) - (b[1]-a[1])*(c[0]-b[0]) ;
/*
 *  Arbitrarily select axis direction so that first non-zero component
 *  or the direction is positive.
 */
sign = 0 ;
if( axis->direction[0] <= 0 )
{
    if( axis->direction[0] < 0 )
         sign = 1 ;
    else if( axis->direction[1] <= 0 )
{
             if( axis->direction[1] < 0 )
                  sign = 1 ;
             else if( axis->direction[2] < 0 )
                      sign = 1 ;
}
}
if( sign )
    for( i = 0 ; i < DIMENSION ; i++ )
        axis->direction[i] = -axis->direction[i] ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += axis->direction[i]*axis->direction[i] ;
r = sqrt(r) ;
for( i = 0 ; i < DIMENSION ; i++ )
    axis->direction[i] /= r ;
if( verbose > 1 ){
    printf( "    axis origin is at (%g,%g,%g)\n", 
        axis->normal[0]*axis->distance, axis->normal[1]*axis->distance, axis->normal[2]*axis->distance ) ;
    printf( "    axis is in the direction (%g,%g,%g)\n", axis->direction[0], axis->direction[1], axis->direction[2] ) ;
    }
return axis ;
}

SYMMETRY_ELEMENT *init_higher_axis( gint ia, gint ib, gint ic )
{
SYMMETRY_ELEMENT * axis ;
gdouble             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
gint                i ;

if( verbose > 0 ) printf( "Trying cn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
StatTotal++ ;
/* Do a quick check of geometry validity */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = Atoms[ia].x[i] - CenterOfSomething[i] ;
    b[i] = Atoms[ib].x[i] - CenterOfSomething[i] ;
    c[i] = Atoms[ic].x[i] - CenterOfSomething[i] ;
    }
if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
    if( verbose > 0 ) printf( "    no coherrent axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_atom ;
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the c%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Improper axes-specific routines.
 *   These are obtained by slight modifications of normal rotation
 *       routines.
 */
void rotate_reflect_atom( SYMMETRY_ELEMENT *axis, OBJECT *from, OBJECT *to )
{
gdouble             x[3], y[3], a[3], b[3], c[3] ;
gdouble             angle = 2*PI/axis->order ;
gdouble             a_sin = sin( angle ) ;
gdouble             a_cos = cos( angle ) ;
gdouble             dot ;
gint                i ;

if (DIMENSION != 3)
  {
  fprintf(stderr, "Catastrophe in rotate_reflect_atom!\n");
  return;
  }

for( i = 0 ; i < 3 ; i++ )
    x[i] = from->x[i] - axis->distance * axis->normal[i] ;
for( i = 0, dot = 0 ; i < 3 ; i++ )
    dot += x[i] * axis->direction[i] ;
for( i = 0 ; i < 3 ; i++ )
    a[i] = axis->direction[i] * dot ;
for( i = 0 ; i < 3 ; i++ )
    b[i] = x[i] - a[i] ;
c[0] = b[1]*axis->direction[2] - b[2]*axis->direction[1] ;
c[1] = b[2]*axis->direction[0] - b[0]*axis->direction[2] ;
c[2] = b[0]*axis->direction[1] - b[1]*axis->direction[0] ;
for( i = 0 ; i < 3 ; i++ )
    y[i] = -a[i] + b[i]*a_cos + c[i]*a_sin ;
for( i = 0 ; i < 3 ; i++ )
    to->x[i] = y[i] + axis->distance * axis->normal[i] ;
to->type = from->type ;
}

SYMMETRY_ELEMENT * init_improper_axis( gint ia, gint ib, gint ic )
{
SYMMETRY_ELEMENT * axis ;
gdouble             a[ DIMENSION ], b[ DIMENSION ], c[ DIMENSION ] ;
gdouble             centerpoint[ DIMENSION ] ;
gdouble             r ;
gint                i ;

if( verbose > 0 ) printf( "Trying sn axis for the triplet (%d,%d,%d)\n", ia, ib, ic ) ;
StatTotal++ ;
/* First, reduce the problem to Cn case */
for( i = 0 ; i < DIMENSION ; i++ ){
    a[i] = Atoms[ia].x[i] - CenterOfSomething[i] ;
    b[i] = Atoms[ib].x[i] - CenterOfSomething[i] ;
    c[i] = Atoms[ic].x[i] - CenterOfSomething[i] ;
    }
for( i = 0, r = 0 ; i < DIMENSION ; i++ ){
    centerpoint[i] = a[i] + c[i] + 2*b[i] ;
    r             += centerpoint[i]*centerpoint[i] ;
    }
r = sqrt(r) ;
if( r <= ToleranceSame ){
    StatEarly++ ;
    if( verbose > 0 ) printf( "    atoms can not define improper axis of the order more than 2\n" ) ;
    return NULL ;
    }
for( i = 0 ; i < DIMENSION ; i++ )
    centerpoint[i] /= r ;
for( i = 0, r = 0 ; i < DIMENSION ; i++ )
    r += centerpoint[i] * b[i] ;
for( i = 0 ; i < DIMENSION ; i++ )
    b[i] = 2*r*centerpoint[i] - b[i] ;
/* Do a quick check of geometry validity */
if( ( axis = init_axis_parameters( a, b, c ) ) == NULL ){
    if( verbose > 0 ) printf( "    no coherrent improper axis is defined by the points\n" ) ;
    return NULL ;
    }
axis->transform_atom = rotate_reflect_atom ;
if( refine_symmetry_element( axis, 1 ) < 0 ){
    if( verbose > 0 ) printf( "    refinement failed for the s%d axis\n", axis->order ) ;
    destroy_symmetry_element( axis ) ;
    return NULL ;
    }
return axis ;
}

/*
 *   Control routines
 */

/****************************************/
/* calculate atom to centroid distances */
/****************************************/
void find_distances(void)
{
gint i;
gdouble r[3];

#if DEBUG
P3VEC("Centroid is ", CenterOfSomething);
#endif

DistanceFromCenter = g_malloc(AtomsCount*sizeof(gdouble));
for (i=0 ; i<AtomsCount ; i++)
  {
  ARR3SET(r, &Atoms[i].x[0]); 
  ARR3SUB(r, CenterOfSomething);
  DistanceFromCenter[i] = VEC3MAGSQ(r);
  } 
}

void find_planes(void)
{
gint i, j;
SYMMETRY_ELEMENT * plane;

plane = init_ultimate_plane() ;
if( plane != NULL ){
    MolecularPlane = plane ;
    PlanesCount++ ;
    Planes = g_realloc(Planes, sizeof(SYMMETRY_ELEMENT *)*PlanesCount);
    Planes[ PlanesCount - 1 ] = plane ;
    }
for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( ( plane = init_mirror_plane( i, j ) ) != NULL ){
            PlanesCount++ ;
            Planes = g_realloc(Planes, sizeof(SYMMETRY_ELEMENT *)*PlanesCount);
            Planes[ PlanesCount - 1 ] = plane ;
            }
        }
    }
}

void find_inversion_centers(void)
{
SYMMETRY_ELEMENT *center;

if ((center = init_inversion_center()) != NULL)
  {
  InversionCenters = g_malloc0(sizeof(SYMMETRY_ELEMENT *));
  InversionCenters[0] = center;
  InversionCentersCount = 1;
  }
}

void find_infinity_axis(void)
{
SYMMETRY_ELEMENT *axis;

if ((axis = init_ultimate_axis()) != NULL)
  {
  NormalAxesCount++;
  NormalAxes = g_realloc(NormalAxes, 
                         sizeof(SYMMETRY_ELEMENT *)*NormalAxesCount);
  NormalAxes[NormalAxesCount-1] = axis;
  }
}

void find_c2_axes(void)
{
gint i, j, k, l, m ;
gdouble center[ DIMENSION ] ;
gdouble *distances;
gdouble r;
SYMMETRY_ELEMENT * axis ;

distances = g_malloc0(AtomsCount*sizeof(gdouble));

for( i = 1 ; i < AtomsCount ; i++ ){
    for( j = 0 ; j < i ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        /*
         *   First, let's try to get it cheap and use CenterOfSomething
         */
        for( k = 0, r = 0 ; k < DIMENSION ; k++ ){
            center[k] = ( Atoms[i].x[k] + Atoms[j].x[k] ) / 2 ;
            r        += pow2( center[k] - CenterOfSomething[k] ) ;
            }
        r = sqrt(r) ;
        if( r > 5*TolerancePrimary ){ /* It's Ok to use CenterOfSomething */
            if( ( axis = init_c2_axis( i, j, CenterOfSomething ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = g_realloc(NormalAxes, sizeof(SYMMETRY_ELEMENT *)
                                       * NormalAxesCount);
                NormalAxes[ NormalAxesCount - 1 ] = axis;
                }
            continue;
            }
        /*
         *  Now, C2 axis can either pass through an atom, or through the
         *  middle of the other pair.
         */
        for( k = 0 ; k < AtomsCount ; k++ ){
            if( ( axis = init_c2_axis( i, j, Atoms[k].x ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = g_realloc(NormalAxes, sizeof(SYMMETRY_ELEMENT *)
                                       * NormalAxesCount);
                NormalAxes[ NormalAxesCount - 1 ] = axis;
                }
            }
        /*
         *  Prepare data for an additional pre-screening check
         */
        for( k = 0 ; k < AtomsCount ; k++ ){
            for( l = 0, r = 0 ; l < DIMENSION ; l++ )
                r += pow2( Atoms[k].x[l] - center[l] ) ;
            distances[k] = sqrt(r) ;
            }
        for( k = 0 ; k < AtomsCount ; k++ ){
            for( l = 0 ; l < AtomsCount ; l++ ){
                if( Atoms[k].type != Atoms[l].type )
                    continue ;
                if( fabs( DistanceFromCenter[k] - DistanceFromCenter[l] ) > TolerancePrimary ||
                    fabs( distances[k] - distances[l] ) > TolerancePrimary )
                        continue ; /* We really need this one to run reasonably fast! */
                for( m = 0 ; m < DIMENSION ; m++ )
                    center[m] = ( Atoms[k].x[m] + Atoms[l].x[m] ) / 2 ;
                if( ( axis = init_c2_axis( i, j, center ) ) != NULL )
                  {
                  NormalAxesCount++ ;
                  NormalAxes = g_realloc(NormalAxes, sizeof(SYMMETRY_ELEMENT * )                                         * NormalAxesCount);
                  NormalAxes[ NormalAxesCount - 1 ] = axis ;
                  }
                }
            }
        }
    }
g_free(distances);
}

void find_higher_axes(void)
{
gint i, j, k ;
SYMMETRY_ELEMENT *axis;

for( i = 0 ; i < AtomsCount ; i++ ){
    for( j = i + 1 ; j < AtomsCount ; j++ ){
        if( Atoms[i].type != Atoms[j].type )
            continue ;
        if( fabs( DistanceFromCenter[i] - DistanceFromCenter[j] ) > TolerancePrimary )
            continue ; /* A very cheap, but quite effective check */
        for( k = 0 ; k < AtomsCount ; k++ ){
            if( Atoms[i].type != Atoms[k].type )
                continue ;
            if( ( fabs( DistanceFromCenter[i] - DistanceFromCenter[k] ) > TolerancePrimary ) ||
                ( fabs( DistanceFromCenter[j] - DistanceFromCenter[k] ) > TolerancePrimary ) )
                    continue ;
            if( ( axis = init_higher_axis( i, j, k ) ) != NULL ){
                NormalAxesCount++ ;
                NormalAxes = g_realloc(NormalAxes, sizeof(SYMMETRY_ELEMENT *)
                                       * NormalAxesCount);
                NormalAxes[ NormalAxesCount - 1 ] = axis ;
                }
            }
        }
    }
}

void find_improper_axes(void)
{
gint                i, j, k ;
SYMMETRY_ELEMENT * axis ;

for( i = 0 ; i < AtomsCount ; i++ ){
  for( j = i + 1 ; j < AtomsCount ; j++ ){
    for( k = 0 ; k < AtomsCount ; k++ ){
      if ((axis = init_improper_axis( i, j, k ) ) != NULL)
        {
        ImproperAxesCount++ ;
        ImproperAxes = g_realloc(ImproperAxes, sizeof(SYMMETRY_ELEMENT *)                                        * ImproperAxesCount);
        ImproperAxes[ ImproperAxesCount - 1 ] = axis ;
        }
      }
    }
  }
}

/*****************/
/* main sequence */
/*****************/
void find_symmetry_elements( void )
{
find_distances();
find_inversion_centers();
find_planes();
find_infinity_axis();
find_c2_axes();
find_higher_axes();
find_improper_axes();
}

gint compare_axes( const void *a, const void *b )
{
SYMMETRY_ELEMENT * axis_a = *(SYMMETRY_ELEMENT**) a ;
SYMMETRY_ELEMENT * axis_b = *(SYMMETRY_ELEMENT**) b ;
gint i, order_a, order_b ;

order_a = axis_a->order ; if( order_a == 0 ) order_a = 10000 ;
order_b = axis_b->order ; if( order_b == 0 ) order_b = 10000 ;
if( ( i = order_b - order_a ) != 0 ) return i ;
if( axis_a->maxdev > axis_b->maxdev ) return -1 ;
if( axis_a->maxdev < axis_b->maxdev ) return  1 ;
return 0 ;
}

void sort_symmetry_elements( void )
{
if( PlanesCount > 1 ){
    qsort( Planes, PlanesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( NormalAxesCount > 1 ){
    qsort( NormalAxes, NormalAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
if( ImproperAxesCount > 1 ){
    qsort( ImproperAxes, ImproperAxesCount, sizeof( SYMMETRY_ELEMENT * ), compare_axes ) ;
    }
}

void summarize_symmetry_elements( void )
{
gint i;

NormalAxesCounts   = g_malloc0((MaxAxisOrder+1)*sizeof(gint));
ImproperAxesCounts = g_malloc0((MaxAxisOrder+1)*sizeof(gint));

for( i = 0 ; i < NormalAxesCount ; i++ )
    NormalAxesCounts[ NormalAxes[i]->order ]++ ;
for( i = 0 ; i < ImproperAxesCount ; i++ )
    ImproperAxesCounts[ ImproperAxes[i]->order ]++ ;
}

/*************************************/
/* update symmetry info in model pak */
/*************************************/
#define  DEBUG_UPDATE_SYMMETRY 0
void update_symmetry(struct model_pak *data)
{
gint i, j , k, n, m, p;
gdouble angle;
GString *code;
struct symop_pak *symops;

/* free the list items */
g_strfreev(data->symmetry.items);

n = PlanesCount + NormalAxesCount + ImproperAxesCount + InversionCentersCount;

#if DEBUG_UPDATE_SYMMETRY
printf("Model has %d symmetry elements:\n", n);
printf("%d Planes,\n", PlanesCount);
printf("%d Rotation axes,\n", NormalAxesCount);
printf("%d Improper rotation axes,\n", ImproperAxesCount);
printf("%d Inversion centers.\n", InversionCentersCount);
#endif

/* nothing found? */
if (!n)
  {
  data->symmetry.num_items = 1;
  data->symmetry.items = g_malloc(2*sizeof(gchar *));
  *(data->symmetry.items) = g_strdup("none");
  *(data->symmetry.items+1) = NULL;
  return;
  }

/* otherwise continue */
g_free((data->symmetry.symops));
data->symmetry.symops = g_malloc(n * sizeof(struct symop_pak));
/* this probably will overallocate (but ensures space for NULL termination) */
data->symmetry.items = g_malloc((n+1) * sizeof(gchar *));

code = g_string_new(NULL);

/* fill out with operations found */
/* j follows number of symop elements/matrices */
symops = data->symmetry.symops;
j = k = 0;
if (InversionCentersCount) 
  {
  g_string_sprintfa(code, "(i) " ) ;
  *((data->symmetry.items)+k++) = g_strdup("Inversion centre");
/* store the symmetry operation */
  (symops+j)->type = INVERSION;
  VEC3SET(&((symops+j)->matrix[0]),-1.0, 0.0, 0.0);
  VEC3SET(&((symops+j)->matrix[3]), 0.0,-1.0, 0.0);
  VEC3SET(&((symops+j)->matrix[6]), 0.0, 0.0,-1.0);
  j++;
  }

/* infinite axes */
if (NormalAxesCounts[0])
  {
  *((data->symmetry.items)+k++) = g_strdup_printf("%d infinite rotation axes",
                                                   NormalAxesCounts[0]);
  g_string_sprintfa(code, "%d*(Cinf) ", NormalAxesCounts[0] ) ;
/*
  j += NormalAxesCounts[0];
*/
  }

p=0;
for (i=MaxAxisOrder ; i>=2 ; i--)
  {
  switch (NormalAxesCounts[i])
    {
    case 0:
      continue;
    case 1:
      g_string_sprintfa(code, "(C%d) ", i );
      break;
    default:
      g_string_sprintfa(code, "%d*(C%d) ", NormalAxesCounts[i], i );
    }
  *((data->symmetry.items)+k++) = g_strdup_printf("%d order %d rotation axes",
                                                    NormalAxesCounts[i],i);
/* store the symmetry operation */
  for (m=0 ; m<NormalAxesCounts[i] ; m++)
    {
    (symops+j)->type = PAXIS;
    angle = 2.0*G_PI / (gdouble) i;

/*
    get_rot_matrix((symops+j)->matrix, PAXIS, &(NormalAxes[p++]->direction[0]), i);
*/
matrix_v_rotation((symops+j)->matrix, &(NormalAxes[p++]->direction[0]), angle);

    j++;
    }
  }

for (i=MaxAxisOrder ; i>=2 ; i--)
  {
/*
  j += ImproperAxesCounts[i];
*/
  switch (ImproperAxesCounts[i])
    {
    case 0:
      continue;
    case 1:
      g_string_sprintfa(code, "(S%d) ", i);
      break;
    default:
      g_string_sprintfa(code, "%d*(S%d) ", ImproperAxesCounts[i], i);
    }
  *((data->symmetry.items)+k++) = g_strdup_printf("%d order %d improper axes",
                                                    ImproperAxesCounts[i],i);
  }

p=0;
switch (PlanesCount) 
  {
  case 0:
    break;
  case 1:
    g_string_sprintfa(code, "(sigma) ");
    *((data->symmetry.items)+k++) = g_strdup_printf("reflection plane");
    break;
  default:
    g_string_sprintfa(code, "%d*(sigma) ", PlanesCount);
    *((data->symmetry.items)+k++) = g_strdup_printf("%d reflection planes",
                                                             PlanesCount);
/* store the symmetry operation */
  for (m=0 ; m<PlanesCount ; m++)
    {
/* NIY */
/* FIXME - infreq. core dump on line below */
/* symops alloc seemed ok, perhaps NormalAxes is not what we want? */
/*
    get_rot_matrix((symops+j)->matrix, PLANE, &(NormalAxes[p++]->direction[0]), i);
    (symops+j)->type = PLANE;
    j++;
*/
    }
  }
#if DEBUG_UPDATE_SYMMETRY
printf("symmetry code: {%s}\n", code->str) ;
printf("Filled out %d items (allocated for %d)\n",j,n);
for (m=0 ; m<j ; m++)
  {
  P3MAT("symop ",(symops+m)->matrix);
  }
#endif
SymmetryCode = g_strdup(code->str);

/* terminate, so that g_strfreev works */
data->symmetry.num_symops = j;
*((data->symmetry.items)+k++) = NULL;
data->symmetry.num_items = k;

g_string_free(code, TRUE);
}

void identify_point_group(struct model_pak *data)
{
gint i;
gint last_matching = -1;
gint matching_count = 0;

for (i=0 ; i<PointGroupsCount ; i++ )
  {
  if (g_ascii_strcasecmp(SymmetryCode, PointGroups[i].symmetry_code) == 0)
    {
    if (PointGroups[i].check() == 1 )
      {
      last_matching = i ;
      matching_count++ ;
      }
#if DEBUG
    else
      printf("It looks very much like %s, but it is not since %s\n", 
              PointGroups[i].group_name, PointGroupRejectionReason);
#endif
    }
  }

#if DEBUG
switch(matching_count)
  {
  case 0:
    printf("No matching point group.\n");
    break;
  case 1:
    printf("Apparent point group: %s\n",PointGroups[last_matching].group_name);
    break;
  default:
    printf("Error in point group lookup.\n");
    return;
  }
#endif

/* update dialog */
if (PointGroups[last_matching].group_name != NULL)
  {
  g_free(data->symmetry.pg_name);
  data->symmetry.pg_name = g_strdup(PointGroups[last_matching].group_name);
  }
else
  {
  g_free(data->symmetry.pg_name);
  data->symmetry.pg_name = g_strdup("unknown");
  }
}

/*
 *  Input/Output
 */
gint get_coords(struct model_pak *data)
{
GSList *list;
struct core_pak *core;

Atoms = g_malloc0(g_slist_length(data->cores) * sizeof(OBJECT));

/* read 'em in */
AtomsCount=0;
for (list=data->cores ; list ; list=g_slist_next(list))
  {
  core = list->data;
  if (core->status & (DELETED | HIDDEN))
    continue;

  Atoms[AtomsCount].type = core->atom_code;
/*
  ARR3SET(Atoms[AtomsCount].x, core->x);
*/
  ARR3SET(Atoms[AtomsCount].x, core->rx);

  AtomsCount++;
  }

/* transfer centroid */
ARR3SET(CenterOfSomething, data->centroid);

return(0);
}

/*************/
/* main call */
/*************/
void gui_symmetry_analyse(void)
{
struct model_pak *data = sysenv.active_model;

/* checks */
if (!data)
  return;
if (!g_slist_length(data->cores))
  return;

/* force all atom view */
unhide_atoms();
data->asym_on = FALSE;

/* init pointers, so we don't free something that wasn't found & alloc'd */
DistanceFromCenter = NULL;
InversionCenters = NULL;
SymmetryCode = NULL;

/* convert formats */
if (get_coords(data))
  {
  printf("Error reading in atomic coordinates\n");
  return;
  }

PlanesCount = NormalAxesCount = ImproperAxesCount = InversionCentersCount = 0;

/* FIXME - globals need init, as successive calls overlap */
find_symmetry_elements();
sort_symmetry_elements();
summarize_symmetry_elements();
#if DEBUG
if (BadOptimization)
  {
  printf("Refinement of some symmetry elements was terminated before"
         "convergence was reached.\n Some symmetry elements may remain"
         "unidentified.\n") ;
  }
#endif
update_symmetry(data);
/* if 0 sym elements, can identify wrongly for some reason */
if (data->symmetry.num_symops)
  identify_point_group(data);

/* free if allocated */
if (DistanceFromCenter)
  g_free(DistanceFromCenter);
if (SymmetryCode)
  g_free(SymmetryCode);
if (InversionCentersCount)
  g_free(InversionCenters);

/* update widget */
gui_active_refresh();
}

/*****************************/
/* hide a symmetry operation */
/*****************************/
#define DEBUG_DELETE_SYMOP 0
gint gui_symmetry_hide(struct model_pak *data, gint num)
{
gint flag, cycle;
gdouble *mat, vec[3], dx[3];
GSList *list1, *list2;
struct core_pak *core1, *core2;

/* get the matrix */
mat = (data->symmetry.symops+num)->matrix;

#if DEBUG_DELETE_SYMOP
P3MAT("",mat);
#endif

for (list1=data->cores ; list1 ; list1=g_slist_next(list1))
  {
  core1 = list1->data;
  if (core1->status & HIDDEN)
    continue;

/* no need for latmat? (ie isolated molecules only) */
  ARR3SET(vec, core1->x);
/* do rotation */
/* multiple applications (until 1st arrive at start again) */
  cycle=0;
  while(!cycle)
    {
    vecmat(mat, vec);

/* which atom of the same type (if any) does it generate? */
    flag=0;
    for (list2=data->cores ; list2 ; list2=g_slist_next(list2))
      {
      core2 = list2->data;
      if (core1->atom_code != core2->atom_code)
        continue;

/* attempt to match rotated i with some j */
      ARR3SET(dx, vec);
      ARR3SUB(dx, core2->x);

      if (VEC3MAGSQ(dx) < POSITION_TOLERANCE)
        {
#if DEBUG_DELETE_SYMOP
printf("[%d : %d]\n",i,j);
#endif
/* we've found at least 1 atom this sym op produces */
        flag++;
/* same core? */
        if (core1 == core2)
          cycle=1;
        else
          core2->status |= HIDDEN;

/* found a match, so we can stop searching */
        break;
        }
      }
/* atom generates nothing -> sym op candidate failed */
    if (!flag)
      return(1);
    }
  }
return(0);
}

/**********************************/
/* reduce to asymetric components */
/**********************************/
void gui_symmetry_toggle(void)
{
gint i;
GSList *list;
struct core_pak *core;
struct model_pak *model = sysenv.active_model;

if (!model)
  return;

model->asym_on ^= 1;

if (model->asym_on)
  {
  switch (model->periodic)
    {
    case 3:
      for (list=model->cores ; list ; list=g_slist_next(list))
        {
        core = list->data;
        if (!core->primary)
          core->status |= HIDDEN;
        }
      break;

    case 0:
      for (i=model->symmetry.num_symops ; i-- ; )
        gui_symmetry_hide(model, i);
      break;

    }
  }
else
  {
/* remove hidden flag */
  for (list=model->cores ; list ; list=g_slist_next(list))
    {
    core = list->data;
    core->status &= ~HIDDEN;
    }
  }

redraw_canvas(SINGLE);
}

/*******************************/
/* symmetry information widget */
/*******************************/
void gui_symmetry_refresh(GtkWidget *box)
{
gint i;
gchar *value;
GtkTreeIter iter;
GtkCellRenderer *renderer;
GtkTreeViewColumn *column;
static GtkListStore *ls=NULL;
static GtkWidget *tv=NULL, *frame, *vbox;
gchar *label[] = {"space group", "system", "atoms",
                      "a", "b", "c", "alpha", "beta", "gamma",
                      "surface area", "volume", NULL};

if (box)
  {
g_assert(ls == NULL);
g_assert(tv == NULL);

/* new tree list */
  ls = gtk_list_store_new(2, G_TYPE_STRING, G_TYPE_STRING);
  tv = gtk_tree_view_new_with_model(GTK_TREE_MODEL(ls));
  gtk_box_pack_start(GTK_BOX(box), tv, TRUE, TRUE, 0);

/* setup cell renderers */
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(" ", renderer, "text", 0, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tv), column);
  renderer = gtk_cell_renderer_text_new();
  column = gtk_tree_view_column_new_with_attributes(" ", renderer, "text", 1, NULL);
  gtk_tree_view_append_column(GTK_TREE_VIEW(tv), column);
  gtk_tree_view_set_headers_visible(GTK_TREE_VIEW(tv), FALSE);

/* currently, allow nothing to be selected */
  gtk_tree_selection_set_mode(gtk_tree_view_get_selection(GTK_TREE_VIEW(tv)),
                              GTK_SELECTION_NONE);

/* actions */
  frame = gtk_frame_new(NULL);
  gtk_box_pack_start(GTK_BOX(box), frame, FALSE, FALSE, 0);
  gtk_container_set_border_width(GTK_CONTAINER(frame), PANEL_SPACING);
  vbox = gtk_vbox_new(TRUE, 1);
  gtk_container_add(GTK_CONTAINER(frame), vbox);

  gui_button_x("Guess Pointgroup ", gui_symmetry_analyse, NULL, vbox);
  gui_button_x("Asymmetric unit ", gui_symmetry_toggle, NULL, vbox);
  }
else
  {
  struct model_pak *model = sysenv.active_model;

g_assert(ls != NULL);
g_assert(tv != NULL);

gtk_list_store_clear(ls);

if (!model)
  return;

/* TODO - put symmetry info in its own pulldown? */
/* Symmetry part */
if (model->periodic)
  {
/* populate with space group info only if 3D */
  if (model->periodic == 3)
    i=0;
  else
    i=3;

  while (label[i])
    {
    value = NULL;
    switch (i)
      {
      case 0:
        value = g_strdup_printf("%s", model->sginfo.spacename);
        break;
      case 1:
        value = g_strdup_printf("%s", model->sginfo.latticename);
        break;
/*
    case 2:
      value = g_strdup_printf("%d (%d)", model->num_atoms, model->num_asym);
      break;
*/
      case 3:
        value = g_strdup_printf("%8.4f", model->pbc[0]);
        break;
      case 4:
        if (model->periodic > 1)
          value = g_strdup_printf("%8.4f", model->pbc[1]);
        break;
      case 5:
        if (model->periodic > 2)
          value = g_strdup_printf("%8.4f", model->pbc[2]);
        break;
      case 6:
        if (model->periodic > 2)
          value = g_strdup_printf("%8.2f", R2D*model->pbc[3]);
        break;
      case 7:
        if (model->periodic > 2)
          value = g_strdup_printf("%8.2f", R2D*model->pbc[4]);
        break;
      case 8:
        if (model->periodic > 1)
          value = g_strdup_printf("%8.2f", R2D*model->pbc[5]);
        break;
      case 9:
        if (model->periodic == 2)
          value = g_strdup_printf("%.4f", model->area);
        break;
      case 10:
        if (model->periodic == 3)
          value = g_strdup_printf("%.2f", model->volume);
        break;
      } 
/* only add if a valid value */
    if (value)
      {
      gtk_list_store_append(ls, &iter);
      gtk_list_store_set(ls, &iter, 0, label[i], -1);
      gtk_list_store_set(ls, &iter, 1, value, -1);
/*
      msb_list = g_slist_prepend(msb_list, value);
*/
      }

    i++;
    }
  }
else
  {
  if (model->id == MORPH)
    {
    i = PG_Index(model->sginfo.pointgroup);

    gtk_list_store_append(ls, &iter);
    gtk_list_store_set(ls, &iter, 0, "point group", 1, PG_Names[i], -1);
    }
  else
    {
/* populate with symmetry info */
    for (i=0 ; i<model->symmetry.num_items ; i++)
      {
      gtk_list_store_append(ls, &iter);
      value = *((model->symmetry.items)+i);
      gtk_list_store_set(ls, &iter, 0, value, -1);
      }
/* point group */
    gtk_list_store_append(ls, &iter);
    gtk_list_store_set(ls, &iter, 0, "point group", 1, model->symmetry.pg_name, -1);
    }
  }

  }

/* fixes data hiding that can occur when columns get resized between models */
gtk_tree_view_columns_autosize(GTK_TREE_VIEW(tv));

}

