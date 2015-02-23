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

gint module_libsiesta_config(gpointer );
void show_siesta_config_page ();
void siesta_dialog();

void file_handler_page_generator(gpointer * , GtkWidget * , struct model_pak *);
void siesta_gui_page_generator(gpointer * , GtkWidget * , struct model_pak *);


void set_basis_set_sz(struct model_pak * );
void set_basis_set_dz(struct model_pak * );
void set_basis_set_szp(struct model_pak * );
void set_basis_set_dzp(struct model_pak * );
void set_basis_set_custom(struct model_pak * );

void set_geom_runtype_sp(struct model_pak * );
void set_geom_runtype_opt(struct model_pak * );
void set_geom_runtype_md(struct model_pak * );
void set_geom_runtype_pc(struct model_pak * );
void set_geom_constantcomp_cp(struct model_pak * );
void set_geom_constantcomp_cv(struct model_pak * );
void kgrid_action (GtkWidget * , gboolean );
void dud_action(GtkWidget *, struct model_pak * );
void set_pulay_sensitive(GtkWidget *, GtkWidget *);
void random_action(GtkWidget *, GtkWidget *);
void optimise_cell_action(GtkWidget *, struct model_pak * );
void long_output_click(GtkWidget *, struct model_pak * );
void zeta_warning (GtkWidget *, struct model_pak * );
void quick_message(GtkWidget *, gchar *);

void set_md_run_type(GtkWidget *, gpointer *);
void set_animate(GtkWidget *, struct model_pak *);


void siesta_prevfreq(GtkWidget *, struct model_pak *);
void siesta_nextfreq(GtkWidget *, struct model_pak *);
    
void siesta_showfreq(GtkWidget *, struct model_pak *);
    
gdouble make_eigvec_freq(gdouble);

gint siesta_timed_animate(struct model_pak *);
void set_animate(GtkWidget *, struct model_pak *);
void siesta_create_timer(struct model_pak *);
void siesta_destroy_timer(struct model_pak *);
void siesta_timeout_cleanup(struct model_pak *);
void siesta_toggle_paused(GtkWidget *, struct model_pak *);
    
void siesta_file_dialog(GtkWidget *, struct model_pak *);
void siesta_file_save_loop(GtkWidget *, struct model_pak *);
    
void siesta_make_runscript(gchar *, gchar *, struct grid_config_pak *);

void siesta_load_animation(GtkWidget *, struct model_pak *);
    
void siesta_animation_dialog_destroy(GtkWidget *, gpointer *);
void siesta_save_n_queue(GtkWidget *, struct model_pak *);


void siesta_show_vect(GtkWidget *, struct model_pak *);
