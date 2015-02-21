
/* fill types for box packing */
enum {FF, TF, FT, TT, LB};

/* generation 2 shortcuts */
void gui_relation_update(gpointer);
void gui_relation_update_widget(gpointer);

GtkWidget *gui_icon_button(const gchar *, const gchar *,
                             gpointer, gpointer,
                             GtkWidget *);
GtkWidget *gui_stock_button(const gchar *,
                              gpointer, gpointer,
                              GtkWidget *);

/* TODO - make auto_check use same protoype as direct_check */
GtkWidget *gui_auto_check(gchar *, gpointer , gpointer , gint *, GtkWidget *);
GtkWidget *gui_direct_check(gchar *, gint *, gpointer , gpointer , GtkWidget *);

/* my shortcut routines for gtk interface construction */
GtkWidget *new_csd(gchar *, gpointer);
GtkWidget *new_check_button(gchar *, gpointer, gpointer, gint, GtkWidget *);

GtkWidget *gui_button(gchar *, gpointer, gpointer, GtkWidget *, gint);
GtkWidget *gui_button_x(gchar *, gpointer, gpointer, GtkWidget *);

void gui_button_label(gchar *, gchar *, gpointer, gpointer, GtkWidget *);

void new_radio_group(gint, GtkWidget *, gint);

void gui_checkbox_refresh(GtkWidget *, GtkWidget *);

GtkWidget *add_radio_button(gchar *, gpointer, gpointer);

GtkWidget *new_spinner(gchar *, gdouble, gdouble, gdouble,
                       gpointer, gpointer, GtkWidget *);

GtkWidget *gui_auto_spin(gchar *, gdouble *, gdouble, gdouble, gdouble,
                           gpointer, gpointer, GtkWidget *);
GtkWidget *gui_direct_spin(gchar *, gdouble *, gdouble, gdouble, gdouble,
                             gpointer, gpointer, GtkWidget *);

GtkWidget *gui_new_spin(gchar *, gdouble, gdouble, gdouble,
                          gpointer, GtkWidget *);

GtkWidget *gui_direct_hscale(gdouble, gdouble, gdouble,
                               gpointer, gpointer, gpointer, GtkWidget *);

GtkWidget *gui_auto_text_label(gchar **);
GtkWidget *gui_auto_int_label(gint *);
GtkWidget *gui_auto_float_label(gdouble *);

GtkWidget *gui_text_window(gchar **, gint);

GtkWidget *gui_text_entry(gchar *, gchar **, gint, gint, GtkWidget *);

GtkWidget *gui_frame_vbox(const gchar *, gint, gint, GtkWidget *);
GtkWidget *gui_frame_hbox(const gchar *, gint, gint, GtkWidget *);

gpointer gui_pulldown_new(const gchar *, GList *, gint, GtkWidget *);
const gchar *gui_pulldown_text(gpointer);

void gui_colour_box(const gchar *, gdouble *, GtkWidget *);


gpointer gui_pd_new(GSList *, gint, gpointer, gpointer);
gchar *gui_pd_text(gpointer);

void gui_hbox_pack(GtkWidget *, gchar *, GtkWidget *, gint);

