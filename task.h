#include <time.h>

/* task modes */
enum {RUNNING, QUEUED, KILLED, COMPLETED, REMOVED};

/************/
/* task pak */
/************/
struct task_pak
{
/* control */
gint pid;
gint ppid;
gint status;
gint parent;
gint child;
gint sister;
gint h_sec;
gint sec;
gint min;
gint hour;
gchar *label;
gchar *message;

gchar *status_file;
FILE *status_fp;
gint status_index;
GString *status_text;
/* JJM DEBUG
#ifdef __WIN32
gpointer start_time;
#else
time_t start_time;
#endif
*/
time_t start_time;

gchar *time;
gdouble pcpu;
gdouble pmem;
/* NEW */
gdouble progress;
gpointer locked_model;

/* main task and arguments (run by the child) */
void (*primary)(gpointer, ...);
gpointer ptr1;
/* cleanup task and arguments (run by the parent) */
void (*cleanup) (gpointer, ...);
gpointer ptr2;
};

/* task control */
gint update_task_info(void);
void task_status_update(struct task_pak *);
void task_dialog(void);
gint task_sync(const gchar *);

void task_queue_init(void);
void task_queue_free(void);

void task_free(gpointer);
void task_new(const gchar *,
              gpointer, gpointer,
              gpointer, gpointer,
              gpointer);

gint exec_gulp(const gchar *, const gchar *);
