WebVis
======

What it is
----------

WebVis is an auxiliary tool providing convenient web interface for the online_vis.
Currently it can:
* display console output and store the log;
* provide auto-completion for commands and file paths;
* send console commands;
* auto-detect filenames and provide downloading URLs;
* auto-convert .dot graphs into .svg when downloading;
* and show then in the window.
What it will be able to do in future:
* provide dialogs for loading environment and other commands;
* dynamically render graphs and update their parts via sigma.js library.

Installing
----------

First, of course, don't forget to build *online_vis* itself – *make rv* in the *assembler* folder.
WebVis is written in Python and Javascript, using Flask_ and jQuery. To launch it, of course you need Python (2.7+ is recommended; not tested for 2.6 and 3.x) and *pip* or *easy_install*. At least one of them is typically pre-installed in both Linux and OS X, so no prerequisites needed. Use it to install the *venv* sandbox:
- sudo easy_install virtualenv
or:
- sudo pip install virtualenv
If you can't *sudo*, ask your administrator.
Then you may run the auto-build script to download and install necessary packages into the sandbox:
- ./build.sh
If you encounter unlikely problems, you may do it manually:
- virtualenv venv
- . venv/bin/activate
- pip install Flask
- pip install Flask-Session
- pip install jsonpickle
- pip install pyparsing

Configuring
-----------

Default configuration is just enough for a successful launch. If you need to change something, all options are set in the *webvis.cfg* file:
- *debug*: when set to *true*, enables debug output in the server stdout. Default is *true*.
- *env*: sets the environment folder, where *online_vis* binaries are located. Default is *../assembler/* (assembler build folder).
- *port*: HTTP server port. Default is 5000.

Running
-------

After installing, you can launch the server:
- ./run.sh
Or, if you dislike scripts:
- . venv/bin/activate #if it wasn't activated already
- python index.py
It listens on *localhost:<port>* (default port is 5000) which can be visited, well, locally. If you plan to install and run it on the remote server and don't want to (or may not) configure ports, you can use an SSH tunnel. Run this command on your machine:
- ssh -L <port>:localhost:<port> <server-login>@<server-ip>
Then you can visit *localhost:<port>* just as if WebVis was installed on your local machine.

Using
-----
There is a small login form at the index page. Enter the unique name to start your session. Names may contain only latin letters and digits. Other logged users are shown to avoid conflicts, which are currently NOT checked. Currently it runs the background *./run rv* process in the *env* folder infinitely. When that process stops working (for example, due to an error), it becomes re-launched.
On the main page, you can send text commands in the prompt as if it was a command line. There is auto-completion for command names and file paths, browsed on the server.
All console output is stored in the session. To stop the session and all related processes, simply click "Log out" on the main page.
When WebVis receives a command response which contains a file or directory path, it automatically converts it into a download hyperlink. If it was a .dot file, a rendered .svg is shown on the page instead. If it is a directory, its content is loaded into the console log, so you can load some of its files. Rendered files are placed in the *static/cache/<username>* folder, so they aren't re-rendered in next time.
You can also traverse the rendered SVG graph by clicking vertices. Then, *draw_vertex* is called and the graph is replaced with a newly rendered one. WARNING: currently that only works for paired graphs and only for the first vertex in the node (due to .dot format limitations).
There is also an experimental feature to render dynamic force-directed graphs which can be expanded with node clicking, but it's hardly useful at present. If you still want to use it, call *draw_dynamic* instead of *draw_vertex*.

Troubleshooting
---------------

*Sending commands, clicking URLS doesn't do anything*
Check if the web server has started. Try to refresh the page and/or re-login. If the server hangs, restart it.

*I get obsolete console log or renders from old sessions*
Clean the cache folder, or run a clean restart:
- ./run.sh clean

*When I send a command, I get a response from the previous one
There is a desynchronization between client and the server. Currently it's only healed with the server restart.

*Command response is messed up, like "No such command `adrw_cnig`"
*Can't launch the server, the port is busy*
There is a background *online_vis* process using the same communication pipes from the previous launch, which was not correctly terminated. Find and kill it with fire^W^W.

*The browser spits out Python error messages*
There is a bug, report that to the developer_.

.. _Flask: http://flask.pocoo.org
.. _developer: mailto:y.s.gorshkov@gmail.com
