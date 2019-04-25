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

#ifdef __WIN32
#define DIR_SEP "\\"
#else
#define DIR_SEP "/"
#endif

#define KEYWORD_LEN 10


#define BOHR_TO_ANGS 0.52917724928
#define HARTREE_TO_EV 27.21162

/* FORTRAN ugliness */
#define READ_RECORD fread(&sysenv.fortran_buffer, sizeof(int), 1, fp)
#define WRITE_RECORD fwrite(&sysenv.fortran_buffer, sizeof(int), 1, fp)

/* useful ADDONs --OVHPA*/
#define __Q(a) #a
#define __SKIP_BLANK(pointer) while(!g_ascii_isgraph(*pointer)&&(*(pointer)!='\0')) pointer++
#define __SKIP_NUM(pointer) while(g_ascii_isdigit(*pointer)&&(*(pointer)!='\0')) pointer++
/*WARN: it is still unsafe to mix fseek/ftell with fgetpos/fsetpos*/
#define __GET_LAST_LINE(fp,buffer) do{\
        fseek(fp,-2,SEEK_END);\
        while(fgetc(fp)!='\n') fseek(fp,-2,SEEK_CUR);\
        fseek(fp,+1,SEEK_CUR);\
        buffer = file_read_line(fp);\
}while(0)

/* enumerated types for POV-Ray colour-style selection */
enum {HSV,	  			/*Colour-wheel style red-green-blue */
	  REDWHITEBLUE		/* Red-fades to white-fades to blue */
};

struct keyword_pak
{
gchar *label;
gint code;
};

extern struct keyword_pak keywords[];

/* main */
void file_init(void);
void file_free(void);
gint fgetline(FILE *, gchar *);
gchar *file_read_line(FILE *);
void file_skip_comment(gint);

gint set_path(const gchar *);
GSList *file_dir_list(const gchar *, gint);
GSList *file_dir_list_pattern(const gchar *, const gchar *);

struct file_pak *get_file_info(gpointer, gint);
void file_load(gchar *, struct model_pak *);
void file_save(void);
gint file_save_as(gchar *, struct model_pak *);
gchar *gun(const gchar *);
gchar *format_value_and_units(gchar *, gint);
gint file_extension_valid(const gchar *);
gchar *file_extension_get(const gchar *);
gint file_byte_size(const gchar *);
gchar *file_find_program(const gchar *);

/* dialog control */
void file_dialog(gchar *, gchar *, gint, 
                 gpointer (gchar *,  struct model_pak *), gint);
void file_load_dialog(void);
void file_save_dialog(void);

/* file writing routines */
gint write_arc(gchar *, struct model_pak *);
gint write_cif(gchar *, struct model_pak *);
gint write_fdf(gchar *, struct model_pak *);
gint write_gulp(gchar *, struct model_pak *);
gint write_gmf(gchar *, struct model_pak *);
gint write_planes(gchar *, struct model_pak *);
gint write_marvin(gchar *, struct model_pak *);
gint write_xml(gchar *, struct model_pak *);
gint write_xtl(gchar *, struct model_pak *);
gint write_xyz(gchar *, struct model_pak *);
gint write_gms(gchar *, struct model_pak *);
gint write_diffax(gchar *, struct model_pak *);
gint write_povray(gchar *, struct model_pak *);
gint write_pdb(gchar *, struct model_pak *);
gint write_gauss(gchar *, struct model_pak *);
gint write_qe(gchar *, struct model_pak *);
gint write_cssr(gchar *, struct model_pak *);
gint write_mol2(gchar *, struct model_pak *);
gint write_dmol(gchar *, struct model_pak *);
gint write_dlpoly(gchar *, struct model_pak *);
gint write_bgf(gchar *, struct model_pak *);
gint write_cgf(gchar *, struct model_pak *);
gint write_dlp(gchar *, struct model_pak *);
gint write_gromacs(gchar *, struct model_pak *);
gint write_castep_cell(gchar *, struct model_pak *);
gint write_meta(gchar *, struct model_pak *);

void write_povray_colour_textures(FILE *, struct model_pak *, int);
void write_sfc_data(FILE *);

gint write_arc_header(FILE *, struct model_pak *);
gint write_arc_frame(FILE *, struct model_pak *);
gint write_trj_header(FILE *, struct model_pak *);
gint write_trj_frame(FILE *, struct model_pak *);

/* file reading routines */
gint read_arc(gchar *, struct model_pak *);
gint read_cif(gchar *, struct model_pak *);
gint read_fdf(gchar *, struct model_pak *);
gint read_xml_vasp(gchar *, struct model_pak *);
gint read_output_uspex(gchar *filename, struct model_pak *model);
gint read_gulp(gchar *, struct model_pak *);
gint read_gulp_output(gchar *, struct model_pak *);
gint read_gmf(gchar *, struct model_pak *);
gint read_planes(gchar *, struct model_pak *);
gint read_marvin(gchar *, struct model_pak *);
gint read_mvnout(gchar *, struct model_pak *);
gint read_sout(gchar *, struct model_pak *);
gint read_xml(gchar *, struct model_pak *);
gint read_xtl(gchar *, struct model_pak *);
gint read_xyz(gchar *, struct model_pak *);
gint read_using_babel(gchar *, struct model_pak *);
gint read_gms(gchar *, struct model_pak *);
gint read_gms_out(gchar *, struct model_pak *);
gint read_diffax(gchar *, struct model_pak *);
gint read_about(gchar *, struct model_pak *);
gint read_nw(gchar *, struct model_pak *);
gint read_nwout(gchar *, struct model_pak *);
gint read_pdb(gchar *, struct model_pak *);
gint read_castep_out(gchar *, struct model_pak *);
gint read_gauss(gchar *, struct model_pak *);
gint read_gauss_out(gchar *, struct model_pak *);
gint read_qe(gchar *, struct model_pak *);
gint read_qe_out(gchar *, struct model_pak *);
gint read_rietica(gchar *, struct model_pak *);
gint read_off(gchar *, struct model_pak *);
gint read_moldy(gchar *, struct model_pak *);
gint read_moldy_restart(gchar *, struct model_pak *);
gint read_cssr(gchar *, struct model_pak *);
gint read_mol2(gchar *, struct model_pak *);
gint read_cel(gchar *, struct model_pak *);
gint read_dmol(gchar *, struct model_pak *);
gint read_dlpoly(gchar *, struct model_pak *);
gint read_bgf(gchar *, struct model_pak *);
gint read_cgf(gchar *, struct model_pak *);
gint read_dlp(gchar *, struct model_pak *);
gint read_gromacs_gro(gchar *, struct model_pak *);

gint project_read(gchar *, struct model_pak *);

gint read_trj_header(FILE *, struct model_pak *);
gint read_trj_frame(FILE *, struct model_pak *, gint);
gint read_arc_frame(FILE *, struct model_pak *);
gint read_sout_frame(FILE *, struct model_pak *);
gint read_xml_vasp_frame(FILE *, struct model_pak *);
gint read_frame_uspex(FILE *vf, struct model_pak *model);
gint read_gms_out_frame(FILE *, struct model_pak *);
gint read_about_frame(FILE *, struct model_pak *);
gint read_nwout_frame(FILE *, struct model_pak *);
gint read_pdb_frame(FILE *, struct model_pak *);
gint read_castep_out_frame(FILE *, struct model_pak *);
gint read_gauss_out_frame(FILE *, struct model_pak *);
gint read_qe_out_frame(FILE *, struct model_pak *);
gint read_xyz_frame(FILE *, struct model_pak *);
gint read_dlpoly_frame(FILE *, struct model_pak *);
gint read_dmol_frame(FILE *, struct model_pak *);

/*NEW: track/update frames (TODO)*/
gint update_frame_uspex(gint idx,struct model_pak *model);

gint load_planes(gchar *, struct model_pak *);

gint mark_trj_frames(struct model_pak *);

void import_planes(gchar *);

void swap_bytes(void *, const gint);

void gdis_blurb(FILE *);

/* parsing */
void capitals(gchar *, gchar *);
gchar **get_tokenized_line(FILE *, gint *);
gint get_keyword_code(const gchar *);

#define find_in_string(a,b) strstr(b,a)
long int fetch_in_file(FILE *vf,const gchar *target);

gint read_frame(FILE *, gint, struct model_pak *);
gint read_raw_frame(FILE *, gint, struct model_pak *);
gint add_frame_offset(FILE *, struct model_pak *);

gint hash_strcmp(gconstpointer, gconstpointer);

GSList *fdf_species_build(struct model_pak *);
gint fdf_species_index(gchar *, GSList *);

GSList *gromacs_read_ff(const gchar *);

