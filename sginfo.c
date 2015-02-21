/*
  Space Group Info's (c) 1994-96 Ralf W. Grosse-Kunstleve
 */

/*
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
#include <ctype.h>
#include <string.h>
#include "gdis.h"


/*
  Macintosh extras (Courtesy Jon Tischler <TischlerJZ@ornl.gov>)
 */
#if defined(__THINK__) || defined(__MWERKS__)
#include <console.h>
#define CONSOLE_LINES   36  /* number of lines to use for console */
#define CONSOLE_COLUMNS 90  /* number of columns to use for console */
#ifdef __MWERKS__
#include <sioux.h>
#endif
#endif


#define AppMalloc(ptr, n) (ptr) = g_malloc((n) * sizeof (*(ptr)))
#define AppFree(ptr, n) free(ptr)


#define SGCOREDEF__
#include "sginfo.h"


#if USE_GS_SI

static int PrimitiveRotMx(const int *CCMx_LP, int *RotMx, const int *CCMx_PL,
                          int deterCCMx_LP)
{
  int       i;
  int       BufMx[9];


  /* Mp = Tlp . Mz . Tpl */

  RotMxMultiply(BufMx, RotMx, CCMx_PL);
  RotMxMultiply(RotMx, CCMx_LP, BufMx);

  for (i = 0; i < 9; i++)
  {
    if (RotMx[i] % deterCCMx_LP) {
      SetSgError("Internal Error: PrimitiveRotMx()");
      return -1;
    }
  }

  for (i = 0; i < 9; i++)
    RotMx[i] /= deterCCMx_LP;

  return 0;
}


static int Find_si(T_SgInfo *SgInfo)
{
  static const int Tab_si_Vector[] =
    {
       1,  0,  0,   0, /*  h      */
       0,  1,  0,   1, /*  k      */
       0,  0,  1,   2, /*  l      */
       1,  1,  0,   0, /*  h+k    */
       1, -1,  0,   0, /*  h-k    */
       0,  1,  1,   1, /*  k+l    */
       0,  1, -1,   1, /*  k-l    */
       1,  0,  1,   1, /*  h+l    */
       1,  0, -1,   1, /*  h-l    */
       1,  1,  1,   0, /*  h+k+l  */
       1,  1, -1,   0, /*  h+k-l  */
       1, -1,  1,   0, /*  h-k+l  */
      -1,  1,  1,   0, /* -h+k+l  */
       2,  1, -1,   0, /*  2h+k-l */
       2, -1,  1,   0, /*  2h-k+l */
      -1,  2,  1,   0, /* -h+2k+l */
       1,  2, -1,   0, /*  h+2k-l */
      -1,  1,  2,   0, /* -h+k+2l */
       1, -1,  2,   0  /*  h-k+2l */
    };

  static int nTab_si_Vector
     = sizeof Tab_si_Vector / sizeof (*Tab_si_Vector) / 4;

  int        deterCCMx_LP, CCMx_PL[9];
  int        i, itabsiv;
  int        nLoopInv, iLoopInv, n_si_v, i_si_v;
  int        n, m, l;
  int        IsFine;
  int        item[3];
  int        R_I[9], si_Buf[9];
  int        iList;
  T_RTMx     *lsmx;
  const int  *tabsiv;


  if (SgInfo->LatticeInfo->Code != 'P')
  {
    deterCCMx_LP = deterRotMx(SgInfo->CCMx_LP);
                 InverseRotMx(SgInfo->CCMx_LP, CCMx_PL);

    if (deterCCMx_LP < 1)
      goto ReturnError;
  }

  nLoopInv = Sg_nLoopInv(SgInfo);

  SgInfo->n_si_Vector = n_si_v = 0;

  for (i = 0; i < 9; i++)
    SgInfo->si_Vector[i] = 0;

  for (i = 0; i < 3; i++)
  {
    SgInfo->si_Modulus[i] = 1;
    item[i] = 1;
  }

  tabsiv = Tab_si_Vector;

  for (itabsiv = 0; itabsiv < nTab_si_Vector; itabsiv++, tabsiv += 4)
  {
    IsFine = 1;
    m = -1;

    for (iList = 0; IsFine && iList < SgInfo->nList; iList++)
    {
      lsmx = &SgInfo->ListSeitzMx[iList];

      for (iLoopInv = 0; IsFine && iLoopInv < nLoopInv; iLoopInv++)
      {
        if (iLoopInv == 0)
          for (i = 0; i < 9; i++)
          {
            if (i % 4) R_I[i] =  lsmx->s.R[i];
            else       R_I[i] =  lsmx->s.R[i] - 1;
          }
        else
          for (i = 0; i < 9; i++)
          {
            if (i % 4) R_I[i] = -lsmx->s.R[i];
            else       R_I[i] = -lsmx->s.R[i] - 1;
          }

        if (SgInfo->LatticeInfo->Code != 'P')
        {
          if (PrimitiveRotMx(SgInfo->CCMx_LP, R_I, CCMx_PL,
                                deterCCMx_LP) < 0)
            return -1;
        }

        for (i = 0; IsFine && i < 3; i++)
        {
          n =  tabsiv[0] * R_I[i * 3 + 0];
          n += tabsiv[1] * R_I[i * 3 + 1];
          n += tabsiv[2] * R_I[i * 3 + 2];
          n = abs(n);

          if (n == 1)
            IsFine = 0;
          else if (m < 2)
            m = n;
          else if (n > 0 && n != m)
            IsFine = 0;
        }
      }
    }

    if (IsFine)
    {
#if DEBUG_Find_si
      fprintf(stdout, "H-Kt %2d %2d %2d   %d\n",
        tabsiv[0], tabsiv[1], tabsiv[2], m);
#endif

      l = tabsiv[3];

      while (item[l] > 1) /* just "if", see break's */
      {
        if (m == item[l]) break;

        if (m == 3 && (   SgInfo->XtalSystem != XS_Trigonal
                       || SgInfo->UniqueDirCode != '=')) break;

        if (m == 4 && (   SgInfo->XtalSystem == XS_Triclinic
                       || SgInfo->XtalSystem == XS_Monoclinic)) break;

        if (m == 2) break;

        /* if (m > 1 || m != 4) break; */

        n_si_v--;
        item[l] = 1;
        break;
      }

      if (item[l] == 1)
      {
        if (itabsiv > 12)
          n_si_v = 0;

        item[l] = m;
        SgInfo->si_Modulus[n_si_v] = m;

        n = n_si_v * 3;
        for (i = 0; i < 3; i++)
          SgInfo->si_Vector[n++] = tabsiv[i];

        n_si_v++;
      }
    }
  }

#if DEBUG_Find_si
  fprintf(stdout, "H-Kt\n");
#endif

  if (SgInfo->LatticeInfo->Code != 'P')
  {
#if DEBUG_Find_si
    for (i = 0; i < n_si_v; i++)
      fprintf(stdout, "H-Kp %2d %2d %2d   %d\n",
        SgInfo->si_Vector[i * 3 + 0],
        SgInfo->si_Vector[i * 3 + 1],
        SgInfo->si_Vector[i * 3 + 2],
        SgInfo->si_Modulus[i]);
    fprintf(stdout, "H-Kp\n");
#endif

    for (i_si_v = 0; i_si_v < n_si_v; i_si_v++)
    {
      for (i = 0; i < 3; i++)
      {
        si_Buf[i_si_v * 3 + i]
          =   SgInfo->si_Vector[i_si_v * 3 + 0] * CCMx_PL[i * 3 + 0]
            + SgInfo->si_Vector[i_si_v * 3 + 1] * CCMx_PL[i * 3 + 1]
            + SgInfo->si_Vector[i_si_v * 3 + 2] * CCMx_PL[i * 3 + 2];
      }
    }

    for (i = 0; i < i_si_v * 3; i++)
    {
      if (si_Buf[i] % deterCCMx_LP)
      {
        n = i / 3; n *= 3;
        fprintf(stdout, " %3d %3d %3d\n",
          si_Buf[n + 0], si_Buf[n + 1], si_Buf[n + 2]);
        goto ReturnError;
      }

      SgInfo->si_Vector[i] = si_Buf[i] / deterCCMx_LP;
    }
  }

  SgInfo->n_si_Vector = n_si_v;
  return n_si_v;

  ReturnError:

  SetSgError("Internal Error: Find_si()");
  return -1;
}


#endif /* USE_GS_SI */

#include <math.h>

typedef struct {
                 double   a, b, c;
                 double   alpha, beta, gamma;
                 double   sa, sb, sg;
                 double   ca, cb, cg;
                 double   v;
                 char     calcs, calcc;
               }
               T_LatticeConstants;


#define PIover180 (PI / 180.0)

#define EpsPI (1.e-6) /* ARBITRARY */


/* ******************************************************************* */


typedef struct
  {
    int                Convention;
    const char         *SgName;
    const T_TabSgName  *InpTSgN;
    const T_TabSgName  *RefTSgN;
    T_RTMx             CBMx, InvCBMx;
  }
  T_SgList;

