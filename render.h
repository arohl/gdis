
/* rendering consts */
#define BALL_SIZE 0.27
#define AMBIENCE_MAX 0.025
#define TRIPLE_SCALE 0.65
#define DOUBLE_SCALE 0.8
#define POVRAY_FONT "timrom.ttf"
/* convert from gdis pixel offset to OpenGL cartesian Angs */
#define PIX2ANG 0.03

/* unified render (povray & openGL) */
enum 
{
STICK, BALL_STICK, CPK, LIQUORICE, POLYHEDRAL, ZONE,
LIGHT_TYPE, DIRECTIONAL, POSITIONAL,
SCALE, HALO_QUALITY, FAST_ROTATION,
ANTIALIAS, SHADOWLESS, DRAW_AXES, PERSPECTIVE, WIRE_FRAME, FOG,
ANIMATE, ANIMATE_TYPE, ANIM_GIF, ANIM_MPEG, ANIM_NAME, ANIM_FRAME,
MORPH_STYLE, MORPH_FINISH
};

/************************/
/* rendering structures */
/************************/

/* TODO - move model to this & replace model->scale with zoom */
struct camera_pak
{
gint mode;         /* free/locked */
gint perspective;  /* true/false */
gdouble fov;       /* perspective field of view */
gdouble zoom;      /* zoom/scale factor */
gdouble x[3];      /* position */
gdouble o[3];      /* orientation vector */
gdouble v[3];      /* viewing vector */
gdouble e[3];      /* crossproduct of viewing and orientation */
gdouble q[4];      /* quaternion modifier */
};

struct light_pak
{
gint type;
gdouble x[3];
gdouble colour[3];
gdouble ambient;
gdouble diffuse;
gdouble specular;
};

/* TODO - replace with general render pak? */
struct povray_pak
{
gint background;   /* colour */
gint camera[3];
gint num_lights;
gint shadowless;
gint animate;
gint atype;
gint axes;
gint delay;
gdouble ambience;
gdouble frad;       /* frame radius */
gchar filename[FILELEN];
gint wire_frame;
gdouble ref_index;
gdouble transmit;
gchar morph_finish[LINELEN];
};

/* prototypes */
void render_make_pipes(GSList **, struct model_pak *);
GSList *render_get_pipes(gint, gint, gint, struct model_pak *);
GSList *render_sort_pipes(GSList *);

void create_polyhedra(struct model_pak *);
void destroy_polyhedra(struct model_pak *);

void povray_task(void);
void povray_exec(gchar *);

