/*
Copyright (C) 2003 by Andrew Lloyd Rohl 

a.rohl@curtin.edu.au

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

/******************************************************************************
 * vector.h
 *
 *	Vector include file containing vector typedef and vector macros
 ******************************************************************************/

/* VECTOR definition */
typedef	struct { gdouble element[3]; } vector;

#define V_E(V,e) ((V).element[e])

#define	V_X(V)	((V).element[0])
#define	V_Y(V)	((V).element[1])
#define	V_Z(V)	((V).element[2])
#define	V_X2(V)	((V).element[0] * (V).element[0])
#define	V_Y2(V)	((V).element[1] * (V).element[1])
#define	V_Z2(V)	((V).element[2] * (V).element[2])

#define V_INIT(A,B,C,D)	{ \
	V_X(A) = B; \
	V_Y(A) = C; \
	V_Z(A) = D; \
	}

#define V_EQUATE(A,B)	{ \
	V_X(A) = V_X(B); \
	V_Y(A) = V_Y(B); \
	V_Z(A) = V_Z(B); \
	}

#define V_SCALER(A,B,C)	{ \
	V_X(A) = (float)(B) * V_X(C); \
	V_Y(A) = (float)(B) * V_Y(C); \
	V_Z(A) = (float)(B) * V_Z(C); \
	}

#define V_ADD(A,b,B,c,C)	{ \
	V_X(A) = (float)b * V_X(B) + (float)c * V_X(C); \
	V_Y(A) = (float)b * V_Y(B) + (float)c * V_Y(C); \
	V_Z(A) = (float)b * V_Z(B) + (float)c * V_Z(C); \
	}

#define V_QADD(A,B,OP,C)	{ \
	V_X(A) = V_X(B) OP V_X(C); \
	V_Y(A) = V_Y(B) OP V_Y(C); \
	V_Z(A) = V_Z(B) OP V_Z(C); \
	}

#define V_0_ASSIGN(A,OP1) {\
	V_X(A) OP1; \
	V_Y(A) OP1; \
	V_Z(A) OP1; \
	}

#define V_1_ASSIGN(A,OP1,B) {\
	V_X(A) OP1 V_X(B); \
	V_Y(A) OP1 V_Y(B); \
	V_Z(A) OP1 V_Z(B); \
	}

#define V_ASSIGN(A,OP1,B,OP2,C)	{ \
	V_X(A) OP1 V_X(B) OP2 V_X(C); \
	V_Y(A) OP1 V_Y(B) OP2 V_Y(C); \
	V_Z(A) OP1 V_Z(B) OP2 V_Z(C); \
	}

#define V_2_ASSIGN(A,OP1,B,OP2,C)	{ \
	V_X(A) OP1 V_X(B) OP2 V_X(C); \
	V_Y(A) OP1 V_Y(B) OP2 V_Y(C); \
	V_Z(A) OP1 V_Z(B) OP2 V_Z(C); \
	}

#define V_3_ASSIGN(A,OP1,B,OP2,C,OP3,D) { \
	V_X(A) OP1 V_X(B) OP2 V_X(C) OP3 V_X(D); \
	V_Y(A) OP1 V_Y(B) OP2 V_Y(C) OP3 V_Y(D); \
	V_Z(A) OP1 V_Z(B) OP2 V_Z(C) OP3 V_Z(D); \
	}
	
#define V_MAGSQ(A) (V_X(A) * V_X(A) + V_Y(A) * V_Y(A) + V_Z(A) * V_Z(A))

#define V_MAG(A) (sqrt(V_MAGSQ(A)))

#define V2_MAGSQ(A) (V_X(A) * V_X(A) + V_Y(A) * V_Y(A))

#define V_DOT(A,B) (V_X(A) * V_X(B) + V_Y(A) * V_Y(B) + V_Z(A) * V_Z(B))

#define V_CROSS(A,B,C) { \
	V_X(A) = V_Y(B) * V_Z(C) - V_Z(B) * V_Y(C); \
	V_Y(A) = V_Z(B) * V_X(C) - V_X(B) * V_Z(C); \
	V_Z(A) = V_X(B) * V_Y(C) - V_Y(B) * V_X(C); \
	}
	
#define V_COSANG(A,B) (V_DOT(A,B) / sqrt(V_MAGSQ(A) * V_MAGSQ(B)))

#define V_ANG(A,B) (acos(V_DOT(A,B) / sqrt(V_MAGSQ(A) * V_MAGSQ(B))))

#define V_SEPARATION(A,B) (sqrt((V_X(A) - V_X(B)) * (V_X(A) - V_X(B)) \
	+ (V_Y(A) - V_Y(B)) * (V_Y(A) - V_Y(B)) \
	+ (V_Z(A) - V_Z(B)) * (V_Z(A) - V_Z(B))))

#define V_SEP_SQ(A,B) ((V_X(A) - V_X(B)) * (V_X(A) - V_X(B)) \
	+ (V_Y(A) - V_Y(B)) * (V_Y(A) - V_Y(B)) \
	+ (V_Z(A) - V_Z(B)) * (V_Z(A) - V_Z(B)))
	
#define V_ZERO(A) (V_X(A) = V_Y(A) = V_Z(A) = 0.0)

