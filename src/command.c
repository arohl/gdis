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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "gdis.h"
#include "file.h"
#include "parse.h"
#include "model.h"
#include "test.h"

/* main data structures */
struct sysenv_pak sysenv;

/***************************/
/* command line processing */
/***************************/
#define DEBUG_COMMAND_MAIN_LOOP 1
void command_main_loop(int argc, char **argv)
{
gint n, status, num_tokens;
gchar *inp, *out, *ext, *tmp, *line, *base, **buff;
gboolean exit=FALSE;
GSList *list, *files, *item, *structures;
struct file_pak *inp_pak, *out_pak;
struct model_pak *model;

if (argc < 2)
  printf("gdis> Commandline interface. Enter ? for help.\n");

/* command line */
do
  {
/* if additional keywords - process as single command line */
  if (argc > 1)
    {
/* convert options into single string for processing */
    buff = g_malloc(argc*sizeof(gchar *));
    for (n=1 ; n<argc ; n++)
      *(buff+n-1) = g_strdup(argv[n]);
    *(buff+argc-1) = NULL;
    line = g_strjoinv(" ", buff);
    g_strfreev(buff);
    exit = TRUE;
    }
  else
    {
/* standard repeating prompt */
    printf("gdis> ");
    line = file_read_line(stdin);
    }

/* quit primitive */
  if (g_ascii_strncasecmp("qu", line, 2) == 0 ||
      g_ascii_strncasecmp("ex", line, 2) == 0)
    exit=TRUE;

/* test primitive */
  if (g_ascii_strncasecmp("te", line, 2) == 0)
    {
    buff = tokenize(line, &num_tokens);
    for (n=1 ; n<num_tokens ; n++)
      test_run(*(buff+n));
    g_strfreev(buff);
    }

/* copy primitive */
  if (g_ascii_strncasecmp("co", line, 2) == 0 || g_ascii_strncasecmp("cp", line, 2) == 0)
    {
    buff = tokenize(line, &num_tokens);

    if (num_tokens > 2)
      {

      out_pak = get_file_info(*(buff+2), BY_EXTENSION);
      if (!out_pak)
        {
        printf("Error: unrecognized output filetype.\n");
        break;
        }
      if (!out_pak->write_file)
        {
        printf("Error: write method for this filetype does not exist.\n");
        break;
        }

/* FIXME! - this is all wrong - should be no dependence on sysenv.cwd for path */
/* since this falls over with something like prompt:>gdis cp /tmp/gok.sot /tmp/dest/gok.meta */
      ext = file_extension_get(*(buff+2));
      files = file_dir_list_pattern(sysenv.cwd, *(buff+1));
      for (list=files ; list ; list=g_slist_next(list))
        {
        inp_pak = get_file_info(list->data, BY_EXTENSION);

        inp = g_build_filename(sysenv.cwd, list->data, NULL);

/* FIXME - changed to correct a problem with Florian's metadata store */
//        base = parse_strip(list->data);
        base = g_strdup(list->data);

        if (inp_pak)
          {
          if (inp_pak->read_file)
            {
            model = model_new();

            status = inp_pak->read_file(inp, model);

            structures = g_slist_find(sysenv.mal, model);

#if DEBUG_COMMAND_MAIN_LOOP
printf("Items loaded = %d [status = %d]\n", g_slist_length(structures), status);
#endif

            if (g_slist_length(structures) > 1)
              {
/* multiple output for input files that contain more than one structure */
              n = 1;
              for (item=structures ; item ; item=g_slist_next(item))
                {
                model = item->data;

                tmp = g_strdup_printf("%s_%d.%s", base, n, ext);
                out = g_build_filename(sysenv.cwd, tmp, NULL);

                status = out_pak->write_file(out, model);

#if DEBUG_COMMAND_MAIN_LOOP
printf("a [%s] -> [%s]\n", inp, out);
#endif
                g_free(out);
                g_free(tmp);
                model_free(model);
                n++;
                }
              }
            else
              {
              tmp = g_strdup_printf("%s.%s", base, ext);

              out = g_build_filename(sysenv.cwd, tmp, NULL);
              g_free(tmp);

              status = out_pak->write_file(out, model);

#if DEBUG_COMMAND_MAIN_LOOP
printf("b [%s] -> [%s]\n", inp, out);
#endif
              g_free(out);
              model_free(model);
              }
            }
          else
            printf("Error: read method for this filetype does not exist.\n");
          }
        else
          {
          printf("Error: unrecognized input filetype.\n");
          }
        g_free(base);
        g_free(inp);
        }


/* CURRENT hack */
if (!files)
  {
/* FIXME */

  }




      g_free(ext);

      free_slist(files);
      }
    else
      printf("Error: insufficient tokens for <copy> command\n");

    g_strfreev(buff);
    }

  if (g_ascii_strncasecmp("ls", line, 2) == 0)
    {
    buff = tokenize(line, &num_tokens);

    if (num_tokens > 1)
      files = file_dir_list_pattern(sysenv.cwd, *(buff+1));
    else
      files = file_dir_list(sysenv.cwd, FALSE);

    for (list=files ; list ; list=g_slist_next(list))
      {
      printf("%s\n", (gchar *) list->data);
      }

    g_strfreev(buff);
    }

  if (g_ascii_strncasecmp("cd", line, 2) == 0)
    {
    buff = tokenize(line, &num_tokens);

    if (num_tokens > 1)
      {
      inp = g_build_path(DIR_SEP, sysenv.cwd, *(buff+1), NULL);
      set_path(inp);
      g_free(inp);
      }

    printf("%s\n", sysenv.cwd);

    g_strfreev(buff);
    }

  if (g_ascii_strncasecmp("pwd", line, 3) == 0)
    {
    printf("%s\n", sysenv.cwd);
    }

  if (g_ascii_strncasecmp("?", line, 1) == 0 || g_ascii_strncasecmp("h", line, 1) == 0)
    {
    printf("GDIS v%4.2f available commands\n\n", VERSION);

    printf("cd <destination>\n");
    printf("    - changes to target directory\n");

    printf("pwd\n");
    printf("    - prints the current working directory\n");

    printf("ls\n");
    printf("    - lists the files in the current working directory\n");

    printf("copy <source> <destination>\n");
    printf("    - convert between filetypes, based on extension (wildcards allowed)\n");
    printf("quit\n");
    printf("    - exit program\n");
    }

/* done -> next line */
  g_free(line);
  }
while (!exit);
}

