WebVis
======

What it is
----------

WebVis is an auxilliary tool providing convenient web interface for the online_vis.
Currently it can:
* display console output ang store the log;
* send console commands;
* auto-detect filenames and provide downloading URLs;
* auto-convert .dot graphs into .png when downloading;
* and show then in the window.
What it will be able to do in future:
* provide dialogs for loading environment and other commands;
* provide auto-completion for commands;
* dynamically render graphs and update their parts via sigma.js library.

Installing
----------

First, of course, don't forget to *make run*.
WebVis is written in Python and Javascript, using Flask_ and jQuery. To launch it, of course you need Python and pip or easy_install. At least one of them is typically pre-installed in both Linux and Mac OS X, so no prerequisities needed.
You may run the auto-build script:
- ./build.sh
If you encounter problems, you may do it manually. The better way is to do it via the *virtualenv* sandbox:
- sudo easy_install virtualenv
or
- sudo pip install virtualenv
Then:
- cd path/to/webvis
- virtualenv venv
- . venv/bin/activate
- pip install Flask
- pip install Flask-Session
- pip install jsonpickle
- pip install pyparsing

Running
-------
After installing, you can launch the server:
- ./run.sh
The same, manual way:
- . venv/bin/activate #if it wasn't activated already
- python index.py
By default (and non-configurable, haha), it listens on localhost:5000.

Using
-----
There is a small login form at the index page. Enter the unique name to start your session (other logged users are shown to avoid conflicts). Currently it runs the background *./run rv* process infinitely. The environment folder, in which binaries must be located, is currently hard-coded to be in the assembler folder (path/to/webvis/../../../) which is a major design flaw and will be fixed someday.
On the main page, you can send text commands in the prompt as if it was a command line. There is only auto-completion for command names.
All console output is stored in the session. To stop the session and all related processes, simply click "Log out" on the main page.
When WebVis receives a command response which contains a file or directory path, it automatically converts it into a download hyperlink. If it was a .dot file, a rendered .svg is shown on the page instead. If it is a directory, its content is loaded into the console log, so you can load some of its files.
You can also traverse the rendered SVG graph by clicking vertices. Then, *draw_vertex* is called and the graph is replaced with a newly rendered one.

Troubleshooting
---------------

*Sending commands, clicking URLS doesn't do anything*
Check if the web server has started. Try to refresh the page and/or re-login.

*Can't launch the server, the port is busy*
Probably a hanging *online_vis* process holds it. Kill it with fire^W^W. Or run the cleaner:
- ./run.sh clean

*The browser spits out Python error messages*
There is a bug, report that to the developer_.

.. _Flask: http://flask.pocoo.org
.. _developer: mailto:y.s.gorshkov@gmail.com
