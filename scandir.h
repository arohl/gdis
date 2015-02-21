/*******************************************************************************
*                                                                              *
*                                   Viewmol                                    *
*                                                                              *
*                              S C A N D I R . C                               *
*                                                                              *
*                 Copyright (c) Joerg-R. Hill, December 2003                   *
*                                                                              *
********************************************************************************
*
* $Id: scandir.h,v 1.1.1.1 2003/08/15 03:56:28 seanfleming Exp $
* $Log: scandir.h,v $
* Revision 1.1.1.1  2003/08/15 03:56:28  seanfleming
* Initial import.
*
* Revision 1.2  2002/02/14 01:48:51  sean
* New configure - should allow non standard include/lib paths to be supplied.
*
* Revision 1.2  2001/09/03 05:02:57  sean
* Marvin's surface generation code added
*
* Revision 1.1.1.1  2001/07/06 02:57:07  andrew
* initial import
*
* Revision 1.2  2003/12/10 15:37:02  jrh
* Release 2.3
*
* Revision 1.1  1999/05/24 01:29:43  jrh
* Initial revision
*
*/

#include <dirent.h>
#include <stdlib.h>
#include <string.h>

/* This function is only required for SunOS, all other supported OS
   have this function in their system library */

#define DEBUG_SCANDIR 1
int scandir(const char *dir, struct dirent ***namelist,
            int (*select)(const struct dirent *),
            int (*compar)(const struct dirent **, const struct dirent **))
{
int j;
  DIR *d;
  struct dirent *entry;
  register int i=0;
  size_t entrysize;

  if ((d=opendir(dir)) == NULL)
     return(-1);

  *namelist=NULL;
  while ((entry=readdir(d)) != NULL)
  {
    if (select == NULL || (select != NULL && (*select)(entry)))
    {
      *namelist=(struct dirent **)realloc((void *)(*namelist),
                 (size_t)((i+1)*sizeof(struct dirent *)));
	if (*namelist == NULL) return(-1);
	entrysize=sizeof(struct dirent)-sizeof(entry->d_name)+strlen(entry->d_name)+1;
	(*namelist)[i]=(struct dirent *)malloc(entrysize);
	if ((*namelist)[i] == NULL) return(-1);

/*
printf("entry: [%s] (%d)\n", entry->d_name, entrysize);
*/

	memcpy((*namelist)[i], entry, entrysize);


	i++;
    }
  }


for (j=0 ; j<10 ; j++)
  printf("<%d>: [%s]\n", j, (*namelist)[j]->d_name);


  if (closedir(d)) return(-1);
  if (i == 0) return(-1);
  if (compar != NULL)
    qsort((void *)(*namelist), (size_t)i, sizeof(struct dirent *), (void *) compar);
    
  return(i);
}

int alphasort(const struct dirent **a, const struct dirent **b)
{
  return(strcmp((*a)->d_name, (*b)->d_name));
}
