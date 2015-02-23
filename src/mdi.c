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

#include "gdis.h"
#include "mdi.h"
#include "mdi_pak.h"

extern struct sysenv_pak sysenv;

/* NEW - solvation core goes here ... */

/* two slightly different functionalities */
/* 1. given solute(s) and solvent - create a mixed box */
/* 2. given a mixed box - insert a molecule (create cavity) */

/* create a new mdi structure */
gpointer mdi_new(void)
{
/* the aim of this is to have references to a number of models, which */
/* can be combined in some fashion to create some solvation model */

return(NULL);
}

void mdi_free(gpointer data)
{
}

/****************************************/
/* create a model from an mdi structure */
/****************************************/
gpointer mdi_model_new(gpointer data)
{
/*
struct mdi_pak *mdi = data;
*/

/* the aim of this is to have references to a number of models, which */
/* can be combined in some fashion to create some solvation model */

/* check models exist */

/* preload coords */

/* combine */

/* initialize new model */

return(NULL);
}

