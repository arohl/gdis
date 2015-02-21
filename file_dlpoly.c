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

/*
 * If you have any questions or comments about DL_POLY file support in GDIS,
 * send it to GDIS mailing list
 * http://lists.sourceforge.net/lists/listinfo/gdis-users
 * (you can also CC to Marcin Wojdyr, search the web for e-mail)
 * 
 * Details about 
 *   The DL_POLY Molecular Simulation Package
 * can be found at http://www.cse.clrc.ac.uk/msi/software/DL_POLY/
 * Following DL_POLY configuration files are (partially) supported:
 *
 *     CONFIG/REVCON                    read-write
 *     HISTORY formatted                read
 *     HISTORY unformatted (binary)     not implemented 
 * HISTORY file provides trajectory and can be very large.
 *
 * To see what features are not supported yet read comments and the code.  
 *
 * When using shell model, there is no way to recognize which "atom"
 * represents core and which shell using only CONFIG/HISTORY file.
 * Convention is used: name of shell "atom" is created with
 * one of following postfixes: -shell _shell -shel _shel -shl _shl -sh _sh
 * eg. Zn-shl or Zn_shel
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "gdis.h"
#include "coords.h"
#include "model.h"
#include "file.h"
#include "parse.h"
#include "scan.h"
#include "matrix.h"
#include "interface.h"

/* main structures */
extern struct sysenv_pak sysenv;
extern struct elem_pak elements[];


/* Replaces all blanks in str with '\0' and put pointers to tokens into tokens.
 * Returns number of tokens (n <= max_split). It doesn't allocates any memory.
 * tokens[] should have length >= max_split. 
 */
int split_string(char *str, int max_split, char **tokens)
{
    char *ptr;
    int counter=0;
    int was_space=1;
    for (ptr=str; *ptr; ++ptr) {
        if (isspace(*ptr)) {
            *ptr=0;
            was_space=1;
        }
        else if (was_space) {
            was_space=0;
            tokens[counter] = ptr;
            counter++;
            if (counter >= max_split)
                break;
        }
    }
    return counter;
}

int get_splitted_line(FILE *fp, int max_split, char **tokens)
{
    gchar line[1024];
    if (fgets(line, 1024, fp))
        return split_string(line, max_split, tokens);
    else
        return 0;
}

void mark_animation_frame(struct model_pak *model, FILE *fp)
{
    fpos_t offset;
    fgetpos(fp, &offset);
    model->frame_list = g_list_append(model->frame_list, 
                                      g_memdup(&offset, sizeof(fpos_t))); 
    model->num_frames++;
}

/**************** file reading ****************/

/* read boundary conditions */
static gint read_dlpoly_pbc(gint imcon, FILE *fp, struct model_pak *model)
{
         /* imcon -- periodic boundary key: 
                imcon       meaning
                0   no periodic boundaries
                1   cubic boundary conditions
                2   orthorhombic boundary conditions
                3   parallelepiped boundary conditions
                4   truncated octahedral boundary conditions
                5   rhombic dodecahedral boundary conditions
                6   x-y parallelogram boundary conditions with
                    no periodicity in the z direction
                7   hexagonal prism boundary conditions
          */
    gint num_tokens, i, j;
    gchar *tokens[64];
    gchar line[1024];
    switch (imcon) {
        case 0:
            model->periodic = 0;
            return 0;
        case 1:
        case 2:
            model->periodic = 3;
            for (i=0; i < 3; ++i) {
                num_tokens = get_splitted_line(fp, 64, tokens);
                if (num_tokens < 3) {
                    gui_text_show(ERROR, "not enought numbers in PBC desc.\n");
                    return 1;
                }
                for (j=0; j<3; ++j)
                    model->latmat[3*i+j] = str_to_float(tokens[j]);
            }
            model->construct_pbc = TRUE;
            matrix_lattice_init(model); 
            return 0;
        case 6:
            /*TODO*/
        case 3:
        case 4:
        case 5:
        case 7:
        default:
            model->periodic = 0;
            /*eat 3 lines*/
            for (i=0; i < 3; ++i) 
                fgets(line, 1024, fp);
            gui_text_show(WARNING, "Sorry, not supported imcon (PBC key). "
                                   "Assuming no periodic boundaries.\n");
            return 0;
    }
}

/* checks if "atom" is a shell in core-shell model */
static gint read_dlpoly_is_shell(gchar *name)
{
    /* check for postfixes: -shell _shell -shel _shel -shl _shl -sh _sh
     * and return postfix length */
    int len = strlen(name);
    if (len > 6 && (g_ascii_strcasecmp(name+len-6, "-shell") == 0
                    || g_ascii_strcasecmp(name+len-6, "_shell") == 0))
        return 6;
    else if (len > 5 && (g_ascii_strcasecmp(name+len-5, "-shel") == 0
                    || g_ascii_strcasecmp(name+len-5, "_shel") == 0))
        return 5;
    else if (len > 4 && (g_ascii_strcasecmp(name+len-4, "-shl") == 0
                    || g_ascii_strcasecmp(name+len-4, "_shl") == 0))
        return 4;
    else if (len > 3 && (g_ascii_strcasecmp(name+len-3, "-sh") == 0
                    || g_ascii_strcasecmp(name+len-3, "_sh") == 0))
        return 3;
    else
        return 0;
}

/* OH -> O, etc */
static void dlpoly_cut_elem_postfix(gchar *elem)
{
    /* don't cut Cl, Ca, Cs, Cr, Co, Cu, Cd, Ce, Os, Hf */
    gchar c = elem[0];
    gchar x = 0;
    if (strlen(elem) <= 1)
        return;
    x = elem[1];
    if (((c=='C' || c=='c') && (x!='l' && x!='L' && x!='a' && x!='A' &&
                                x!='s' && x!='S' && x!='r' && x!='R' &&
                                x!='o' && x!='O' && x!='u' && x!='U' &&
                                x!='d' && x!='D' && x!='e' && x!='E'))
        || ((c=='O' || c=='o') && (x!='s' && x!='S')) 
        || ((c=='H' || c=='h') && (x!='f' && x!='F')))
        elem[1] = '\0';
}

static gint read_dlpoly_atoms(gint levcfg, FILE *fp, struct model_pak *model)
{
    /* Possible values of levcfg: 1 - only coordinates in file
                                  2 - coordinates and velocities
                                  3 - coordinates, velocities and forces 
    */
    gchar elem[9];
    gchar *tokens[64];
    gint rec=0, num_tokens, sl=0;
    gdouble *x=0, *v=0;
    struct core_pak *core=NULL;
    struct shel_pak *shell=NULL;
    GSList *core_list, *shell_list;

    model->cores = g_slist_reverse(model->cores);
    model->shels = g_slist_reverse(model->shels);
    core_list = model->cores;
    shell_list = model->shels;

    while ((num_tokens = get_splitted_line(fp, 64, tokens))) {
        if (g_ascii_strncasecmp(tokens[0], "timestep", 8) == 0) {
            break;
        }
        switch (rec) {
            case 0: /* record i */
                if (num_tokens < 1) {
                    gui_text_show(ERROR, "Unexpected syntax in DL_POLY file"
                                         " (in record i).\n");
                    return 1;
                }
                sl = read_dlpoly_is_shell(tokens[0]);
                if (sl) {
                    g_strlcpy(elem, tokens[0], 9);
                    elem[strlen(elem) - sl] = '\0';
                    dlpoly_cut_elem_postfix(elem);
                    if (shell_list) {
                        shell = (struct shel_pak *) shell_list->data;
                        shell_list = g_slist_next(shell_list);
                    }
                    else {
                        shell = shell_new(elem, tokens[0], model);
                        model->shels = g_slist_prepend(model->shels, shell);
                    }
                    x = shell->x;
                    v = shell->v;
                }
                else { //core
                    g_strlcpy(elem, tokens[0], 9);
                    dlpoly_cut_elem_postfix(elem);
                    if (core_list) {
                        core = (struct core_pak *) core_list->data;
                        core_list = g_slist_next(core_list);
                    }
                    else {
                        core = core_new(elem, tokens[0], model);
                        model->cores = g_slist_prepend(model->cores, core);
                    }
                    x = core->x;
                    v = core->v;
                }
                break;
            case 1: /* record ii */
                if (num_tokens < 3) {
                    gui_text_show(ERROR, "Unexpected syntax in DL_POLY file"
                                         " (in record ii).\n");
                    return 1;
                }
                VEC4SET(x, str_to_float(tokens[0]), 
                           str_to_float(tokens[1]), 
                           str_to_float(tokens[2]), 
                           1.0);
                break;
            case 2: /* record iii */
                if (num_tokens < 3) {
                    gui_text_show(ERROR, "Unexpected syntax in DL_POLY file"
                                         " (in record iii).\n");
                    return 1;
                }
                VEC3SET(v, str_to_float(tokens[0]), 
                           str_to_float(tokens[1]), 
                           str_to_float(tokens[2]));
                break;
            case 3: /* record iv */
                if (num_tokens < 3) {
                    gui_text_show(ERROR, "Unexpected syntax in DL_POLY file"
                                         " (in record iv).\n");
                    return 1;
                }
                /* no force info in GDIS */
                break;
            default: /* never comes here */
                printf("ERROR\n");
        }
        rec = (rec+1) % (levcfg+2); /* rec = 0, 1, ..., levcfg+1, 0, 1, ... */
    }
    model->cores = g_slist_reverse(model->cores);
    model->shels = g_slist_reverse(model->shels);
    printf("# cores: %i; # shells: %i\n", g_slist_length(model->cores),
                                          g_slist_length(model->shels));
    return 0;
}


gint read_dlpoly(gchar *filename, struct model_pak *model)
{
    gint num_tokens, history;
    gchar line[1024], *tokens[64];
    FILE *fp;
    gint levcfg, /* tells if velocities and forces are in file*/
         imcon;  /* periodic boundary key */
    int i;

    fp = fopen(filename, "rt");
    if (!fp)
        return 1;

    /* first line is a title */
    if (!fgets(line, 1024, fp))
        return 2;
    /* cut trailing blanks */
    for (i=strlen(line)-1; i>0 && isspace(line[i]); --i)
        line[i] = 0;
    model->title = g_strdup(line);

    /* second line */
    num_tokens = get_splitted_line(fp, 64, tokens);
    if (num_tokens < 2) 
        return 3;
    levcfg = (gint) str_to_float(tokens[0]);
    imcon = (gint) str_to_float(tokens[1]);

    /* 3rd line -- is this CONFIG or HISTORY? */
    if (!fgets(line, 1024, fp))
        return 4;
    history = (g_ascii_strncasecmp(line, "timestep", 8) == 0);

    if (history) { /* HISTORY */
        /* read first frame */
        num_tokens = split_string(line, 64, tokens);
        if (num_tokens >= 6)
            ; /* TODO read current nstep and time-step */
        read_dlpoly_pbc(imcon, fp, model);
        read_dlpoly_atoms(levcfg, fp, model);

        /* mark all frames offsets */
        model->num_frames = 0;
        fseek(fp, 0, SEEK_SET);
        /* TODO optimize it */
        while (fgets(line, 1024, fp)) {
            if (g_ascii_strncasecmp(line, "timestep", 8) == 0) {
                fseek(fp, -strlen(line), SEEK_CUR);
                mark_animation_frame(model, fp);
                fseek(fp, strlen(line), SEEK_CUR);
            }
        }
    }
    else { /* CONFIG/REVCON */
        if (num_tokens >= 4)
            ; /*TODO read current nstep and time-step */
        read_dlpoly_pbc(imcon, fp, model);
        read_dlpoly_atoms(levcfg, fp, model);
    }

    /* model setup */
    strcpy(model->filename, filename);
    g_free(model->basename);
    model->basename = parse_strip(filename);
    model->fractional = FALSE;
    model_prep(model);

    return 0;
}

/******** animation -- overwrite frame ********/
gint read_dlpoly_frame(FILE *fp, struct model_pak *model)
{
    gint num_tokens;
    gint imcon, levcfg;
    gchar *tokens[64];

    puts("in read_dlpoly_frame()");
    g_assert(fp != NULL);

    num_tokens = get_splitted_line(fp, 64, tokens);
    if (num_tokens < 6)
        return 11;
    g_assert(g_ascii_strncasecmp(tokens[0], "timestep", 8) == 0);

    levcfg = (gint) str_to_float(tokens[3]);
    imcon = (gint) str_to_float(tokens[4]);
    read_dlpoly_pbc(imcon, fp, model);
    read_dlpoly_atoms(levcfg, fp, model);
    return 0;
}

/**************** file writing ****************/
gint write_dlpoly(gchar *filename, struct model_pak *model)
{
    GSList *list;
    struct core_pak *core;
    struct shel_pak *shell;
    FILE *fp;
    gint levcfg=0, imcon=0;
    gdouble x[3];
    int i;
    
    /* checks */
    g_return_val_if_fail(model != NULL, 1);
    g_return_val_if_fail(filename != NULL, 2);
    
    /* open the file */
    fp = fopen(filename,"wt");
    if (!fp)
      return(3);
    
    /* first line is a title */
    fprintf(fp, "%s\n", model->title);
    
    if (model->periodic == 0) {
        imcon = 0; //no PBC
    }
    else if (model->periodic == 3) {
        if (model->latmat[0*3+1] == 0.  && model->latmat[0*3+2] == 0.
                && model->latmat[1*3+0] == 0. && model->latmat[1*3+2] == 0.
                && model->latmat[2*3+0] == 0. && model->latmat[2*3+1] == 0.) {
            if (model->latmat[0*3+0] == model->latmat[1*3+1]
                    && model->latmat[0*3+0] == model->latmat[2*3+2])
                imcon = 1; //cubic
            else
                imcon = 2; //orthorhombic
        }
        /*TODO other cases */
    }
    /* TODO when levcfg should be 1? */
    fprintf(fp, "%10i%10i\n", levcfg, imcon);
    /* TODO write current nstep and time-step */
    if (imcon > 0) { /*write PBC*/
        for (i=0; i<3; ++i)
            fprintf(fp, "%20.10f%20.10f%20.10f\n", model->latmat[3*i], 
                                                   model->latmat[3*i+1],
                                                   model->latmat[3*i+2]);
    }
    for (list=model->cores ; list ; list=g_slist_next(list)) {
        core = list->data;
        if (core->status & DELETED)
            continue;

        fprintf(fp, "%-8s\n", core->atom_label); 
    
        ARR3SET(x, core->x);
        vecmat(model->latmat, x);
        fprintf(fp, "%20.8f%20.8f%20.8f\n", x[0], x[1], x[2]);
        if (levcfg > 0)
            fprintf(fp,"%20.8f%20.8f%20.8f\n",core->v[0],core->v[1],core->v[2]);
    }
    /*TODO order of atoms (shells and cores) should be preserved,
     * now if there are shells, they go at the end */
    for (list=model->shels ; list ; list=g_slist_next(list)) {
        shell = list->data;
        if (shell->status & DELETED)
            continue;

        fprintf(fp, "%-8s\n", shell->shell_label); 
    
        ARR3SET(x, shell->x);
        vecmat(model->latmat, x);
        fprintf(fp, "%20.8f%20.8f%20.8f\n", x[0], x[1], x[2]);
        if (levcfg > 0)
            fprintf(fp, "%20.8f%20.8f%20.8f\n", shell->v[0], shell->v[1], 
                                                                shell->v[2]);
    }
    fclose(fp);
    return 0;
}

