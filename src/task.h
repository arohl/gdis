#include <time.h>

/* task modes */
enum {RUNNING, QUEUED, KILLED, COMPLETED, REMOVED};

/************/
/* task pak */
/************/
struct task_pak
{
/* control */
gboolean is_async;/*we NEED async tasks to get correct task PID*/
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
long int status_fp_pos;
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
gint task_async(const gchar *command,pid_t *pid);
gint task_sync_now(const gchar *);

void task_queue_init(void);
void task_queue_free(void);

void task_free(gpointer);
void task_new(const gchar *,
              gpointer, gpointer,
              gpointer, gpointer,
              gpointer);
void task_system_new(const gchar *label,
              gpointer func1, gpointer arg1,
              gpointer func2, gpointer arg2,
              gpointer model);
gint exec_gulp(const gchar *, const gchar *);
gint bg_exec_gulp(const gchar *, const gchar *, struct task_pak *);


