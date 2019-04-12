# Getting started
# We start this tutorial with a very simple example that creats an empty window of size 400x200 pixels and adds a button to it
using Gtk

win = GtkWindow("My First Gtk.jl Program", 400, 200)

b = GtkButton("Click Me")
push!(win,b)

showall(win)

# We will now extend the example to let the button actually do something.
using Gtk

win = GtkWindow("My First Gtk.jl Program", 400, 200)

b = GtkButton("Click Me")
push!(win,b)

function on_button_clicked(w)
  println("The button has been clicked")
end
signal_connect(on_button_clicked, b, "clicked")

showall(win)

# Properties
# If you're following along, you probably noticed that creating win caused quite a lot of output:
GtkWindowLeaf(name="", parent, width-request=-1, height-request=-1, visible=TRUE, sensitive=TRUE, app-paintable=FALSE, can-focus=FALSE, has-focus=FALSE, is-focus=FALSE, can-default=FALSE, has-default=FALSE, receives-default=FALSE, composite-child=FALSE, style, events=0, no-show-all=FALSE, has-tooltip=FALSE, tooltip-markup=NULL, tooltip-text=NULL, window, double-buffered=TRUE, halign=GTK_ALIGN_FILL, valign=GTK_ALIGN_FILL, margin-left=0, margin-right=0, margin-top=0, margin-bottom=0, margin=0, hexpand=FALSE, vexpand=FALSE, hexpand-set=FALSE, vexpand-set=FALSE, expand=FALSE, border-width=0, resize-mode=GTK_RESIZE_QUEUE, child, type=GTK_WINDOW_TOPLEVEL, title="My window", role=NULL, resizable=TRUE, modal=FALSE, window-position=GTK_WIN_POS_NONE, default-width=-1, default-height=-1, destroy-with-parent=FALSE, hide-titlebar-when-maximized=FALSE, icon, icon-name=NULL, screen, type-hint=GDK_WINDOW_TYPE_HINT_NORMAL, skip-taskbar-hint=FALSE, skip-pager-hint=FALSE, urgency-hint=FALSE, accept-focus=TRUE, focus-on-map=TRUE, decorated=TRUE, deletable=TRUE, gravity=GDK_GRAVITY_NORTH_WEST, transient-for, attached-to, opacity=1.000000, has-resize-grip=TRUE, resize-grip-visible=TRUE, application, ubuntu-no-proxy=FALSE, is-active=FALSE, has-toplevel-focus=FALSE, startup-id, mnemonics-visible=TRUE, focus-visible=TRUE, )

# Layout
# The most simple layout widget is the GtkBox. It can be either be horizontally or vertical aligned. It allow to add an arbitrary number of widgets.
win = GtkWindow("New title")
hbox = GtkBox(:h)  # :h makes a horizontal layout, :v a vertical layout
push!(win, hbox)
cancel = GtkButton("Cancel")
ok = GtkButton("OK")
set_gtk_property!(hbox,:expand,ok,true)
set_gtk_property!(hbox,:spacing,10)
push!(hbox, cancel)
push!(hbox, ok)
showall(win)
#More generally, you can arrange items in a grid:
win = GtkWindow("A new window")
g = GtkGrid()
a = GtkEntry()  # a widget for entering text
set_gtk_property!(a, :text, "This is Gtk!")
b = GtkCheckButton("Check me!")
c = GtkScale(false, 0:10)     # a slider

# Now let's place these graphical elements into the Grid:
g[1,1] = a    # Cartesian coordinates, g[x,y]
g[2,1] = b
g[1:2,2] = c  # spans both columns
set_gtk_property!(g, :column_homogeneous, true)
set_gtk_property!(g, :column_spacing, 15)  # introduce a 15-pixel gap between columns
push!(win, g)
showall(win)

# Signals and Callbacks
# Let's try a simple example:
b = GtkButton("Press me")
win = GtkWindow(b, "Callbacks")
showall(win)

function button_clicked_callback(widget)
    println(widget, " was clicked!")
end

id = signal_connect(button_clicked_callback, b, "clicked")
# Using Julia's do syntax, the exact same code could alternatively be written as
b = GtkButton("Press me")
win = GtkWindow(b, "Callbacks")
id = signal_connect(b, "clicked") do widget
     println(widget, " was clicked!")
end
signal_handler_disconnect(b, id)
id2 = signal_connect(b, "clicked") do widget
    println("\"", get_gtk_property(widget,:label,String), "\" was clicked!")
end
showall(win)
