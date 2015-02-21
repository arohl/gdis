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
#include <stdlib.h>

#include "gdis.h"
#include "grid.h"
#include "parse.h"
#ifdef WITH_GRISU
#include "grisu_client.h"
#endif

/************************************************/
/* print list consisting of pointers to strings */
/************************************************/
void test_print_text_list(GSList *list)
{
GSList *item;

for (item=list ; item ; item=g_slist_next(item))
  {
  printf("[%s]\n", (gchar *) item->data);
  }
}

/*******************************/
/* run tests on grisu WS calls */
/*******************************/
void test_grisu(void)
{
#ifdef WITH_GRISU
gint passed=0, total=0;
GSList *list;
gchar *password, *name;
gchar *contents;
gsize length;
gchar *username, buff[999];
gchar *dir;

printf("testing grisu ...\n");

grisu_job_names();

//printf("\nTEST - EXPERIMENTAL\n");

//grisu_df();
//grisu_file_list("gsiftp://ngdata.ivec.org/store/ng03/grid-admin/C_AU_O_APACGrid_OU_iVEC_CN_Sean_Fleming/");

/* NEW - uploading & downloading works */
/* TODO - stress testing / binary content */
//grisu_file_download("gsiftp://ngdata.ivec.org/store/ng03/grid-admin/C_AU_O_APACGrid_OU_iVEC_CN_Sean_Fleming/gibb_best3.got", "/home/sean/prog/gdis/gotcha.got");


//grisu_file_upload("/home/sean/prog/gdis/gdis_manual.txt", "gsiftp://ngdata.ivec.org/store/ng03/grid-admin/C_AU_O_APACGrid_OU_iVEC_CN_Sean_Fleming/upload.txt");
//grisu_file_download("gsiftp://ngdata.ivec.org/store/ng03/grid-admin/C_AU_O_APACGrid_OU_iVEC_CN_Sean_Fleming/upload.txt", "TODO_use_this");

/*
printf("URL for directory listing: ");
if (fgets(buff, sizeof(buff), stdin) != NULL)
  dir = parse_strip_newline(buff);
else
  {
  printf("Input error.\n");
  return;
  }
grisu_file_list(dir);
*/

/* NEW - job submission now works!!! */
/*
if (g_file_get_contents("test.xml", &contents, &length, NULL))
  {
grisu_auth_set(TRUE);
name = grisu_job_submit(contents);
  }
else
  {
printf("Failed to open sample job description XML file in cwd.\n");
name = NULL;
  }

if (name)
  printf("submitted [%s] : success!\n", name);
else
  printf("failed.\n");

grisu_ps();
*/

return;


grisu_auth_set(FALSE);
total++;
printf("\nTEST 1 - Submit location query (no auth)\n");
list = grisu_submit_find("mrbayes");
if (list)
  {
  test_print_text_list(list);
  passed++;
  }

total++;
printf("\nTEST 2 - Submit location query (with auth)\n");
grisu_auth_set(TRUE);
list = grisu_submit_find("mrbayes");
if (list)
  {
  test_print_text_list(list);
  passed++;
  }

total++;
printf("\nTEST 3 - Retrieve available VOs\n");
grisu_auth_set(TRUE);
list = grisu_fqans_get();
if (list)
  {
  test_print_text_list(list);
  passed++;
  }

/* TODO - upload/download file */
/*
grisu_auth_set(TRUE);
name = grisu_job_submit("rubbish");
*/

/*
grisu_auth_set(TRUE);
grisu_absolute_job_dir("gok", "ng2.ivec.org", "/ARCS/Startup");
*/

/* TODO - submit job */
/*
grisu_upload("home/sean/prog/gdis/README.grisu", "grisu-jobs-dir/gok/test.txt");
*/

printf("grisu tests passed: %d/%d\n", passed, total);

#else
printf("grisu component not compiled.\n");
#endif
}

void test_run(gchar *name)
{

#ifdef WITH_GRISU
if (g_ascii_strncasecmp(name, "gri", 3) == 0)
  {
/* DEBUG - re-use existing (already uploaded) proxy */
/*
printf("Enter your MYPROXY credential username: ");
if (fgets(buff, sizeof(buff), stdin) != NULL)
  {
  username = parse_strip_newline(buff);
  }
else
  {
  printf("Input error.\n");
  return;
  }
printf("Enter your MYPROXY credential password: ");
password = parse_getline_hidden();
grisu_keypair_set(username, password);
*/

grisu_keypair_set("someguy", "qetuo135");
grid_setup();

/* cmd line prompt for passwd ... */
/*
printf("Please enter your GRID passphrase: ");
password = parse_getline_hidden();
grisu_keypair_set(NULL, password);
grid_connect(GRID_LOCALPROXY);
g_free(password);
*/

  test_grisu();
  grisu_stop();
  }
#endif

}


