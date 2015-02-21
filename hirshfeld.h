/*
 *  hirshfeld.h
 *  gdis
 *
 *  Created by Joshua McKinnon on 24/08/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

/* Hirshfeld surface prototypes */
void hfs_init(void);
gdouble hfs_calc_wf(gdouble, gdouble, gdouble, struct model_pak *, GSList *);
GSList *hfs_calc_env(GSList *, struct model_pak *);
int hfs_calc_normals(GSList *, struct model_pak *, GSList *);
int ssatoms_calc_normals(GSList *points, struct model_pak *model);
gdouble ssatoms_calc_density(gdouble, gdouble, gdouble, struct model_pak *);

int hfs_calulation_limits(GSList *selection, struct model_pak *, gdouble *min, gdouble *max);
