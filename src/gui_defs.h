/*
Copyright (C) 2018 by Okadome Valencia

hubert.valencia _at_ imass.nagoya-u.ac.jp

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

/* graphic define interface */
/**************************************************** 
 * By using compiler defines,  this is aimed to be an
 * easy to upgrade interface design: by just changing
 * the content here will upgrade all files gui_* that
 * use a specific GTK interface...            --OVHPA
 ****************************************************/

/* DEFINES: interface */
#define GUI_SPACE 0
#define GUI_OBJ GtkWidget
#define GUI_NOTE GtkNotebook
#define GUI_TEXTVIEW GtkTextView
#define GUI_TEXTVIEW_BUFFER GtkTextBuffer

/*****************/
/* BOX INTERFACE */
/*****************/
/*Create a new frame (frame) in box (box)*/
/*Q: add? gtk_container_set_border_width(GTK_CONTAINER(frame), GUI_SPACE);*/
#define GUI_FRAME_BOX(box,frame) do{\
	frame = gtk_frame_new(NULL);\
	gtk_box_pack_start(GTK_BOX(box),frame,FALSE,FALSE,0);\
}while(0)
/*Set a new horizontal box (hbox) in a box (box).*/
#define GUI_LINE_BOX(box,hbox) do{\
        hbox = gtk_hbox_new(FALSE, 0);\
        gtk_container_set_border_width(GTK_CONTAINER(hbox), GUI_SPACE);\
        gtk_box_pack_start(GTK_BOX(box), hbox, FALSE, FALSE, 0);\
}while(0)
/*Set a new label (label) with text ("text") in a box (box).*/
#define GUI_LABEL_BOX(box,label,text) do{\
        label = gtk_label_new(text);\
        gtk_box_pack_start(GTK_BOX(box), label, FALSE, FALSE, 0);\
}while(0)
/*Set a new text entry (entry) with initial value ("value") expanding (wilde=TRUE) or not (wilde=FALSE) in an horizontal box (hbox).*/
#define GUI_TEXT_ENTRY(hbox,entry,value,wilde) do{\
	gchar *_text;\
        entry=gtk_entry_new();\
	_text=g_strdup_printf("%s",value);\
        if(value!=NULL) gtk_entry_set_text(GTK_ENTRY(entry),_text);\
        gtk_box_pack_start(GTK_BOX(hbox), entry, wilde, wilde, 0);\
	g_free(_text);\
}while(0)
/*Set a new entry (entry) with initial value (value) of format ("format") with a width of (size) characters, in an horizontal box (hbox).*/
#define GUI_ENTRY(hbox,entry,value,format,size) do{\
	gchar *_text;\
        entry=gtk_entry_new();\
	_text=g_strdup_printf(format,value);\
        gtk_entry_set_text(GTK_ENTRY(entry),_text);\
        gtk_entry_set_width_chars (GTK_ENTRY(entry),size);\
        gtk_box_pack_start(GTK_BOX(hbox), entry, FALSE, FALSE, 0);\
	g_free(_text);\
}while(0)
/*Set a new separator (separator) in an horizontal box (hbox).*/
#define GUI_NEW_SEPARATOR(hbox,separator) do{\
        separator=gtk_hseparator_new();\
        gtk_box_pack_start(GTK_BOX(hbox), separator, TRUE, FALSE, 0);\
}while(0)
/*Set a new combo (combo), with an uneditable link to GLIST (list) which is FREED after link, with an initial value ("default_text"), and connected to function (function), in an horizontal box (hbox).*/
/* THIS GTK COMBO USE IS DEPRECATED! */
#define GUI_REG_COMBO(hbox,combo,list,default_text,function) do{\
        combo = gtk_combo_new();\
_Pragma ("GCC warning \"use of GTK COMBO interface is deprecated!\"");\
        gtk_entry_set_editable(GTK_ENTRY(GTK_COMBO(combo)->entry), FALSE);\
        gtk_combo_set_popdown_strings(GTK_COMBO(combo), list);\
        g_list_free(list);\
        gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);\
        gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(combo)->entry),default_text);\
        g_signal_connect(GTK_OBJECT(GTK_COMBO(combo)->entry), "changed",GTK_SIGNAL_FUNC(function),data);\
}while(0)
/*Create can open button (button) connected on click to function (function) in box (box).*/
#define GUI_OPEN_BUTTON_BOX(box,button,function) do{\
	button=gtk_button_new_from_stock(GTK_STOCK_OPEN);\
	gtk_box_pack_end(GTK_BOX(box), button, FALSE, FALSE, 0);\
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(function),NULL);\
}while(0)
/*Create an gtk-absurdly complex uneditable, scrollable textview (tv) with textbuffer (tv_buffer), in a box (box).*/
#define GUI_TEXTVIEW_BOX(box,tv,tv_buffer) do{\
	GtkWidget *_absurd = gtk_scrolled_window_new(NULL, NULL);\
        tv = gtk_text_view_new();\
        tv_buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(tv));\
        gtk_text_view_set_editable(GTK_TEXT_VIEW(tv),FALSE);\
	gtk_widget_set_sensitive(tv,FALSE);\
        gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(_absurd),GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);\
        gtk_container_add(GTK_CONTAINER(_absurd),tv);\
        gtk_container_set_border_width(GTK_CONTAINER(_absurd),1);\
        gtk_box_pack_start(GTK_BOX(box),_absurd,TRUE,TRUE,0);\
}while(0)
/*Set text ("text") in textbuffer (tv_buffer)*/
#define GUI_TEXTVIEW_SET(tv_buffer,text) do{\
	gtk_text_buffer_set_text(tv_buffer,text,-1);\
}while(0)
/*Insert text ("text") in textbuffer (tv_buffer) - at cursor.*/
#define GUI_TEXTVIEW_INSERT(tv_buffer,text) do{\
	gtk_text_buffer_insert_at_cursor(tv_buffer,text,-1);\
}while(0)
/**********************/
/* NOTEBOOK INTERFACE */
/**********************/
/*create a new page (page) with label ("caption") in notebook (notebook).*/
#define GUI_PAGE_NOTE(notebook,page,caption) do{\
	GtkWidget *_label = gtk_label_new(caption);\
	page = gtk_vbox_new(FALSE, GUI_SPACE);\
	gtk_notebook_append_page(GTK_NOTEBOOK(notebook),page,_label);\
}while(0)
/*Create a new frame (frame) with label ("caption") in a notebook page (page).*/
#define GUI_FRAME_NOTE(page,frame,caption) do{\
	frame = gtk_frame_new(caption);\
	gtk_box_pack_start(GTK_BOX(page),frame,FALSE,FALSE,0);\
}while(0)
/*create a new table (table) in a notebook page (page).*/
#define GUI_TABLE_NOTE(page,table,row,col) do{\
	table=gtk_table_new(row, col, FALSE);\
	gtk_box_pack_start(GTK_BOX(page),table,TRUE,TRUE,0);\
}while(0)
/*Connect a switch-page even of a notebook (notebook) to a function (function) passing data (data).*/
#define GUI_PAGE_CHANGE(notebook,function,data) do{\
	g_signal_connect(GTK_NOTEBOOK(notebook),"switch-page",GTK_SIGNAL_FUNC(function),data);\
}while(0)
/*return page number*/
#define GUI_NOTE_PAGE_NUMBER(notebook,number) do{\
	number = gtk_notebook_get_current_page(GTK_NOTEBOOK(notebook));\
}while(0)
/*******************/
/* FRAME INTERFACE */
/*******************/
/*Create a new vertical box (vbox) in a frame (frame).*/
#define GUI_VBOX_FRAME(frame,vbox) do{\
        vbox = gtk_vbox_new(FALSE, GUI_SPACE);\
        gtk_container_add(GTK_CONTAINER(frame), vbox);\
}while(0)
/*Create a new unexpanded (row) x (col) table (table) in a frame (frame).*/
#define GUI_TABLE_FRAME(frame,table,row,col) do{\
	table = gtk_table_new(row, col, FALSE);\
	gtk_container_add(GTK_CONTAINER(frame), table);\
}while(0)
/*Create a new toptab notebook (notebook) with border in a frame (frame).*/
#define GUI_NOTE_FRAME(frame,notebook) do{\
	notebook = gtk_notebook_new();\
	gtk_notebook_set_tab_pos(GTK_NOTEBOOK(notebook), GTK_POS_TOP);\
	gtk_container_add(GTK_CONTAINER(frame), notebook);\
	gtk_notebook_set_show_border(GTK_NOTEBOOK(notebook), TRUE);\
}while(0)
/************************************************************************************/
/* TABLE INTERFACE:                                                                 */
/* 	all of the *_TABLE() defines set a new table cell @l,r,t,b in table (table).*/
/************************************************************************************/
#define _ATTACH_TABLE(table,widget,l,r,t,b) gtk_table_attach_defaults(GTK_TABLE(table),widget,l,r,t,b)
//gtk_table_attach(GTK_TABLE(table),widget,l,r,t,b,GTK_FILL,GTK_FILL,0,0)
/*Create a new cell with:
 * 	 a label with text ("text")*/
#define GUI_LABEL_TABLE(table,text,l,r,t,b) do{\
	GtkWidget *_label = gtk_label_new(text);\
	_ATTACH_TABLE(table,_label,l,r,t,b);\
}while(0)
/*Create a new frame with the text ("text")*/
#define GUI_FRAME_TABLE(table,text,l,r,t,b) do{\
	GtkWidget *_frame = gtk_frame_new(text);\
	_ATTACH_TABLE(table,_frame,l,r,t,b);\
}while(0)
/*Create a new cell with:
 *	 a boxed label (created with GUI_LABEL_BOX) with text ("caption"),
 *	 an expanded text entry (entry, created with GUI_TEXT_ENTRY) with a default value ("value").*/
#define GUI_TEXT_TABLE(table,entry,value,caption,l,r,t,b) do{\
	GtkWidget *_label;\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
        GUI_LABEL_BOX(_hbox,_label,caption);\
	GUI_TEXT_ENTRY(_hbox,entry,value,TRUE);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a boxed label (created with GUI_LABEL_BOX) with text ("caption"),
 * 	 a boxed separator,
 * 	 a sized-8 boxed entry (entry, created with GUI_ENTRY) with an initial value (value) of format ("format").*/
#define GUI_ENTRY_TABLE(table,entry,value,format,caption,l,r,t,b) do{\
	GtkWidget *_separator;\
	GtkWidget *_label;\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	GUI_LABEL_BOX(_hbox,_label,caption);\
	GUI_NEW_SEPARATOR(_hbox,_separator);\
	GUI_ENTRY(_hbox,entry,value,format,8);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a new check button (check) with an initial value (value) connected to function (function) with text ("caption").*/
#define GUI_CHECK_TABLE(table,check,value,function,caption,l,r,t,b) do{\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	check = gui_direct_check(caption,&(value),function,NULL,_hbox);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a boxed label (created with GUI_LABEL_BOX) with text ("caption"),
 * 	 a boxed separator,
 * 	 a boxed combo (combo,created with GUI_REG_COMBO) linked to GLIST (list) which is FREED after GUI_REG_COMBO.
 * 	 */
/* THIS GTK COMBO USE IS DEPRECATED! */
#define GUI_COMBO_TABLE(table,combo,default_text,function,caption,l,r,t,b) do{\
	GtkWidget *_separator;\
	GtkWidget *_label;\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
_Pragma ("GCC warning \"use of GTK COMBO interface is deprecated!\"");\
	GUI_LABEL_BOX(_hbox,_label,caption);\
	GUI_NEW_SEPARATOR(_hbox,_separator);\
	GUI_REG_COMBO(hbox,combo,list,default_text,function);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a separator (not boxed!).
 * a.k.a. an empty cell.*/
#define GUI_SEPARATOR_TABLE(l,r,t,b) do{\
	GtkWidget *_separator = gtk_hseparator_new();\
	gtk_table_attach_defaults(GTK_TABLE(table),_separator,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a boxed label with text ("caption"),
 * 	 an extended boxed combobox which is set as uneditable.*/
/* THIS IS THE PROPER (but not simplier) REPLACEMENT FOR GTK COMBO DEPRECATED INTERFACE! */
#define GUI_COMBOBOX_TABLE(table,combobox,caption,l,r,t,b) do{\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	GtkWidget *_label = gtk_label_new(caption);\
        gtk_box_pack_start(GTK_BOX(_hbox),_label,FALSE,FALSE,0);\
        combobox=gtk_combo_box_text_new_with_entry();\
        gtk_entry_set_editable(GTK_ENTRY(gtk_bin_get_child(GTK_BIN(combobox))),FALSE);\
        gtk_box_pack_start(GTK_BOX(_hbox),combobox,TRUE,TRUE,0);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Add a new element with text ("text") in the combobox (combobox).
 * This have to be called AFTER GUI_COMBOBOX.*/
#define GUI_COMBOBOX_ADD(combobox,text) do{\
	gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox),text);\
}while(0)
/*Set a default value (default_value) of the combobox (combobox), and register the function (function) for a changed value in combobox.
 * This have to be called AFTER GUI_COMBOBOX and preferentially after GUI_COMBOBOX_ADD added some elements to combobox.
 * If GUI_COMBOBOX_ADD is called after, it may trigger the function (function) each time.*/
#define GUI_COMBOBOX_SETUP(combobox,default_value,function) do{\
	gtk_combo_box_set_active(GTK_COMBO_BOX(combobox),default_value);\
	if((function)!=NULL) g_signal_connect(GTK_COMBO_BOX_TEXT(combobox),"changed",GTK_SIGNAL_FUNC(function),data);\
}while(0)
/*Only set a default value (default_value) for the combobox (combobox).
 * Can be call anytime, and is equivalent to GUI_COMBOBOX_SETUP(combobox,default_value,NULL).*/
#define GUI_COMBOBOX_SET(combobox,default_value) do{\
	gtk_combo_box_set_active(GTK_COMBO_BOX(combobox),default_value);\
}while(0)
/*Get combobox (combobox) active index (index).*/
#define GUI_COMBOBOX_GET(combobox,index) do{\
	index=gtk_combo_box_get_active(GTK_COMBO_BOX(combobox));\
}while(0)
/*Get combobox (combobox) active text ("text")*/
#define GUI_COMBOBOX_GET_TEXT(combobox,text) do{\
	text=gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(combobox));\
}while(0)
/*Insert combobox (combobox) text ("text") at position index (index).*/
/* NEW: use non-deprecated gtk_combo_box_text_insert_text */
#define GUI_COMBOBOX_ADD_TEXT(combobox,index,text) do{\
	gtk_combo_box_text_insert_text(GTK_COMBO_BOX_TEXT(combobox),index,text);\
}while(0)
/*Prepend text ("text") in combobox (combobox).*/
#define GUI_COMBOBOX_PREPEND_TEXT(combobox,text) do{\
	gtk_combo_box_text_prepend_text(GTK_COMBO_BOX_TEXT(combobox),text);\
}while(0)
/*Delete combobox (combobox) text at position index (index)*/
#define GUI_COMBOBOX_DEL(combobox,index) do{\
	gtk_combo_box_text_remove(GTK_COMBO_BOX_TEXT(combobox),index);\
}while(0)
/*Wipe the combobox (combobox)*/
#define GUI_COMBOBOX_WIPE(combobox) do{\
	gtk_list_store_clear(GTK_LIST_STORE(gtk_combo_box_get_model(GTK_COMBO_BOX(combobox))));\
}while(0)
/*Create a new cell with:
 * 	 a boxed spin button, which default value is 1 and interval is [1,100], and labeled with text ("caption").
 * Due to limitation of GTK, value has to be of double type and not integer (event though it make little sense).*/
#define GUI_SPIN_TABLE(table,spin,value,function,caption,l,r,t,b) do{\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	spin = gui_direct_spin(caption,&(value),1.0,100.0,1.0,function,NULL,_hbox);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Set spin (spin) within range [min,max] (min,max).*/
#define GUI_SPIN_RANGE(spin,min,max) do{\
	gtk_spin_button_set_range(GTK_SPIN_BUTTON(spin),min,max);\
}while(0)
/*Set spin (spin) at value (value).*/
#define GUI_SPIN_SET(spin,value) do{\
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin),value);\
}while(0)
/*Create a new cell with:
 * 	 an open button (button) connected on click to function (function).*/
#define GUI_OPEN_BUTTON_TABLE(table,button,function,l,r,t,b) do{\
	button=gtk_button_new_from_stock(GTK_STOCK_OPEN);\
	_ATTACH_TABLE(table,button,l,r,t,b);\
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(function),NULL);\
}while(0)
/*Create a new cell with:
 * 	 an apply button (button) connected on click to function (function).*/
#define GUI_APPLY_BUTTON_TABLE(table,button,function,l,r,t,b) do{\
	button=gtk_button_new_from_stock(GTK_STOCK_APPLY);\
	_ATTACH_TABLE(table,button,l,r,t,b);\
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(function),NULL);\
}while(0)
/*Create a new cell with:
 * 	 a delete button (button) connected on click to function (function).*/
#define GUI_DELETE_BUTTON_TABLE(table,button,function,l,r,t,b) do{\
	button=gtk_button_new_from_stock(GTK_STOCK_DELETE);\
	_ATTACH_TABLE(table,button,l,r,t,b);\
	g_signal_connect(GTK_OBJECT(button),"clicked",GTK_SIGNAL_FUNC(function),NULL);\
}while(0)
/*Create a new cell with:
 * 	 a boxed button of type APPLY connected on click to a function (function_1),
 * 	 a boxed button of type DELETE connected on click to a function (function_2).*/
#define GUI_2BUTTONS_TABLE(table,function_1,function_2,l,r,t,b) do{\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	GtkWidget *_sep;\
	GUI_NEW_SEPARATOR(_hbox,_sep);\
	gui_stock_button(GTK_STOCK_APPLY,function_1,NULL,_hbox);\
	gui_stock_button(GTK_STOCK_DELETE,function_2,NULL,_hbox);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*left aligned labels*/
/*Create a new cell with:
 * 	 a boxed label with text ("caption"),
 * 	 a separator.
 * a.k.a. a left aligned label.*/
#define GUI_LEFT_LABEL_TABLE(table,caption,l,r,t,b) do{\
	GtkWidget *_separator;\
	GtkWidget *_hbox = gtk_hbox_new(FALSE, 0);\
	GtkWidget *_label = gtk_label_new(caption);\
	gtk_box_pack_start(GTK_BOX(_hbox),_label,FALSE,FALSE,0);\
	GUI_NEW_SEPARATOR(_hbox,_separator);\
	_ATTACH_TABLE(table,_hbox,l,r,t,b);\
}while(0)
/*Create a new cell with:
 * 	 a radio button (radio_1) with caption ("caption_1") connected to function (pointer_1)
 * 	 a radio button (radio_2) with caption ("caption_2") connected to function (pointer_2)*/
#define GUI_2RADIO_TABLE(table,radio_1,caption_1,pointer_1,radio_2,caption_2,pointer_2,l,r,t,b) do{\
	GtkWidget *_vbox = gtk_vbox_new(FALSE,0);\
	new_radio_group(0,_vbox,FF);\
	radio_1 = add_radio_button(caption_1,(gpointer)pointer_1,NULL);\
	radio_2 = add_radio_button(caption_2,(gpointer)pointer_2,NULL);\
	_ATTACH_TABLE(table,_vbox,l,r,t,b);\
}while(0)
/********************/
/* GENERAL COMMANDS */
/********************/
/*create a tooltip with text ("text") for the widget (widget).*/
#define GUI_TOOLTIP(widget,text) gtk_widget_set_tooltip_text(widget,text)
/*registering value of entry (entry) into variable (value), according to format (format).*/
#define GUI_REG_VAL(entry,value,format) do{\
        sscanf(gtk_entry_get_text(GTK_ENTRY(entry)),format,&(value));\
}while(0)
/*set the entry (entry) text ("text")*/
#define GUI_ENTRY_TEXT(entry,text) gtk_entry_set_text(GTK_ENTRY(entry),text);
/*connect entry (entry) on change with the function (function) passing data (data).*/
#define GUI_ENTRY_CHANGE(entry,function,data) g_signal_connect(GTK_OBJECT(GTK_ENTRY(entry)),"changed",GTK_SIGNAL_FUNC(function),data)
/*get the text ("text") of an entry (entry).*/
/* NEW: GTK2 says that return text is const and shall not be freed, hence the g_strdup
 * so text HAVE TO be freed after use.*/
#define GUI_ENTRY_GET_TEXT(entry,text) do{\
	text=g_strdup(gtk_entry_get_text(GTK_ENTRY(entry)));\
}while(0)
/*get the lenght of text in entry (entry)*/
#define GUI_ENTRY_LENGHT(entry) (gtk_entry_get_text_length(GTK_ENTRY(entry)))
/*set sensitivity of a widget*/
#define GUI_LOCK(widget) gtk_widget_set_sensitive(widget,FALSE)
#define GUI_UNLOCK(widget) gtk_widget_set_sensitive(widget,TRUE)
/*set a toggle button (button) ON and OFF*/
#define GUI_TOGGLE_ON(button) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),TRUE)
#define GUI_TOGGLE_OFF(button) gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button),FALSE)
/*insert a frame (frame) in a dialog (window) using an implicit vertical box.*/
#define GUI_FRAME_WINDOW(window,frame) do{\
	GUI_FRAME_BOX(GTK_DIALOG(window)->vbox,frame);\
}while(0)
/*add a save action button (button) connected to function (function) passing data (data) on the dialog (window).*/
#define GUI_SAVE_ACTION(window,button,function,data) do{\
	button=gui_stock_button(GTK_STOCK_SAVE,function,data,GTK_DIALOG(window)->action_area);\
}while(0)
/*add an execute action button (button) connected to function (function) passing data (data) on the dialog (window).*/ 
#define GUI_EXEC_ACTION(window,button,function,data) do{\
	button=gui_stock_button(GTK_STOCK_EXECUTE,function,data,GTK_DIALOG(window)->action_area);\
}while(0)
/*add a close action button (button) connected to function (function) passing data (data) on the dialog (window).*/
#define GUI_CLOSE_ACTION(window,button,function,data) do{\
	button=gui_stock_button(GTK_STOCK_CLOSE,function,data,GTK_DIALOG(window)->action_area);\
}while(0)
/*prepare an open-dialog (open_dialog) in the window (window) with the title (title) and an implicit filter with pattern (filter_pattern) and filter title (filter_caption).*/
#define GUI_PREPARE_OPEN_DIALOG(window,open_dialog,title,filter_pattern,filter_caption) do{\
	GtkFileFilter *_filter = gtk_file_filter_new();\
	gtk_file_filter_add_pattern(_filter,filter_pattern);\
	if(filter_caption!=NULL) gtk_file_filter_set_name (_filter,filter_caption);\
	open_dialog = gtk_file_chooser_dialog_new(title,GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);\
	gtk_file_chooser_add_filter (GTK_FILE_CHOOSER(open_dialog),_filter);\
}while(0)
/*prepare an open-folder dialog (open_dialog) in the window (window) with the title (title).*/
#define GUI_PREPARE_OPEN_FOLDER(window,open_dialog,title) do{\
	open_dialog = gtk_file_chooser_dialog_new(title,GTK_WINDOW(window),GTK_FILE_CHOOSER_ACTION_OPEN,GTK_STOCK_CANCEL,GTK_RESPONSE_CANCEL,GTK_STOCK_OPEN,GTK_RESPONSE_ACCEPT,NULL);\
	gtk_file_chooser_set_action(GTK_FILE_CHOOSER(open_dialog),GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER);\
}while(0)
/*run the previously open-dialog (open_dialog) and either:
 * 	i) get selected filename (filename) and signal having an answer (have_answer);
 * 	OR
 * 	ii) set filename to NULL and signal absence of answer (have_answer). */
#define GUI_OPEN_DIALOG_RUN(open_dialog,have_answer,filename) do{\
	if(gtk_dialog_run(GTK_DIALOG(open_dialog)) == GTK_RESPONSE_ACCEPT) {\
		have_answer=TRUE;\
		filename=gtk_file_chooser_get_filename (GTK_FILE_CHOOSER (file_chooser));\
	}else{\
		have_answer=FALSE;\
		filename=NULL;\
	}\
}while(0)
/*destroy the open-dialog (open_dialog)*/
#define GUI_KILL_OPEN_DIALOG(open_dialog) do{\
	gtk_widget_destroy (GTK_WIDGET(open_dialog));\
}while(0)
/*display the dialog (window).*/
#define GUI_SHOW(window) do{\
	gtk_widget_show_all(window);\
}while(0)
/*close window (window)*/
#define GUI_CLOSE(window) do{\
	gtk_widget_destroy(window);\
}while(0)
