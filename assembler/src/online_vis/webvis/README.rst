WebVis
======

What it is
----------

WebVis is an auxilliary tool providing convenient web interface for the online_vis.
Currently it can:
* display console output;
* send console commands;
* auto-detect filenames and provide downloading URLs;
* auto-convert .dot graphs into .png when downloading.
What it will be able to do in future:
* store logs of commands and files;
* provide dialogs for loading environment and other commands;
* show rendered graphs in the box;
* dynamically render graphs via d3js library.

Installing
----------

WebVis is written in Python and Javascript, using Flask_ and jQuery . To launch it, of course you need Python and pip or easy_install. They are typically pre-installed in both Linux and Mac OS X, so no prerequisities needed. Then you need to install Flask via the VirtualEnv:
- sudo easy_install virtualenv
or
- sudo pip virtualenv
Then:
- cd path/to/webvis
- virtualenv venv
- . venv/bin/activate
- pip install Flask
- pip install Flask-Session

Running
-------
After installing, you can launch the server:
- . venv/bin/activate #if it wasn't activated already
- python index.py
By default (and non-configurable, haha), it listens on localhost:5000.

Using
-----
When the index page is visited, a global session is launched automatically. Currently it runs the background *./run rv* process infinitely. The environment folder, in which binaries must be located, is currently hard-coded to be in the assembler folder (path/to/webvis/../../../).
All console output is stored in the session. To stop the session and all related processes, simply click "Log out" on the main page.
On the main page, you can send text commands in the prompt as if it was a command line. No auto-completion, though.
When WebVis receives a command response which contains a filepath, it automatically converts it into a download hyperlink. If it was a .dot file, a rendered .png is downloaded instead. Probably it will be a configurable option.

.. _Flask: http://flask.pocoo.org
