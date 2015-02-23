
enum {IDENTITY, INVERSION, PLANE, PAXIS, IAXIS, REFLECTION, ALIGNMENT, LATMAT};

/* dim 2 */
#define ARR2SET(vec, a) {*vec = *a ; *(vec+1) = *(a+1);}
#define ARR2ADD(vec, a) {*vec += *a ; *(vec+1) += *(a+1);}
#define VEC2ADD(vec, a, b) {*vec += a ; *(vec+1) += b;}
#define VEC2SET(vec, a, b) {*vec = a ; *(vec+1) = b;}

/* dim 3 */
#define VEC3MUL(vec,m) {*vec *= m ; *(vec+1) *= m ; *(vec+2) *= m;}
#define VEC2MUL(vec,m) {*vec *= m ; *(vec+1) *= m;}
#define VEC3MAGSQ(vec) (*vec * *vec + *(vec+1) * *(vec+1) + \
                        *(vec+2) * *(vec+2))
#define VEC3MAG(vec) sqrt(VEC3MAGSQ(vec))
#define VEC3ABS(a) {*a*=fabs(*a); *(a+1)*=fabs(*(a+1)); *(a+2)*=fabs(*(a+2));}

#define VEC3SET(a,x,y,z) {*a = x; *(a+1) = y; *(a+2) = z;}
#define VEC3ADD(a,x,y,z) {*a += x; *(a+1) += y; *(a+2) += z;}
#define VEC3SUB(a,x,y,z) {*a -= x; *(a+1) -= y; *(a+2) -= z;}
#define ARR3SET(a,b) {*a = *b; *(a+1) = *(b+1); *(a+2) = *(b+2);}
#define ARR3ADD(a,b) {*a += *b; *(a+1) += *(b+1); *(a+2) += *(b+2);}
#define ARR3SUB(a,b) {*a -= *b; *(a+1) -= *(b+1); *(a+2) -= *(b+2);}
#define ARR3MUL(a,b) {*a *= *b; *(a+1) *= *(b+1); *(a+2) *= *(b+2);}

#define DOT3PROD(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

/* dim 4 */
#define VEC4SET(a,x,y,z,t) {*a = x; *(a+1) = y; *(a+2) = z; *(a+3) = t;}
#define VEC4MUL(v,m) {*v *= m; *(v+1) *= m; *(v+2) *= m; *(v+3) *= m;}
#define ARR4SET(a,b) {*a=*b; *(a+1)=*(b+1); *(a+2)=*(b+2); *(a+3)=*(b+3);}

#define P4MAT(s,m) {printf("%s\n|%f, %f, %f, %f|\n|%f, %f, %f, %f|\n\
|%f, %f, %f, %f|\n|%f, %f, %f, %f|\n",s,*m,*(m+1),*(m+2),*(m+3),*(m+4),\
*(m+5),*(m+6),*(m+7),*(m+8),*(m+9),*(m+10),*(m+11),*(m+12),*(m+13),*(m+14),\
*(m+15));}


/* matrix/vector manipulation */
void vecmat(gdouble *, gdouble *);
void vectmat(gdouble *, gdouble *);
void vec4mat(gdouble *, gdouble *);
void matmat(gdouble *, gdouble *);
void mat4mat(gdouble *, gdouble *);
gint matrix_invert(gdouble *);
void matrix_transpose(gdouble *);
void matrix_transpose_44(gdouble *);

void calc_norm(gdouble *, gdouble *, gdouble *, gdouble *);

void matrix_lattice_init(struct model_pak *);
void matrix_lattice_new(gdouble *, struct model_pak *);
void matrix_identity(gdouble *);
void matrix_rotation(gdouble *, gdouble , gint);

void matrix_z_alignment(gdouble *, gdouble *);
void matrix_v_alignment(gdouble *, gdouble *, gdouble *);
void matrix_v_rotation(gdouble *, gdouble *, gdouble);
void matrix_v_reflection(gdouble *, gdouble *);

void matrix_relative_rotation(gdouble *, gdouble, gint, struct model_pak *);

void crossprod(gdouble *, gdouble *, gdouble *);
void proj_vop(gdouble *, gdouble *, gdouble *);

gdouble calc_sep(gdouble *, gdouble *);
gdouble calc_volume(gdouble *);
gint normalize(gdouble *, gint);
gdouble magnitude(gdouble *, gint);
gdouble via(gdouble *, gdouble *, gint);

gdouble matrix_determinant(gdouble *);

gint matrix_is_inversion(gdouble *);
gint matrix_is_identity(gdouble *);
gint matrix_is_z_rotation(gdouble *);
gint matrix_order(gdouble *);

