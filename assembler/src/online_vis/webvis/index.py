#Flask imports
import flask
from flask import Flask, session, request
from flask.ext.session import Session
#Aux imports
import re
from urllib import quote, unquote
#System imports
from os import listdir
from os.path import isfile, join
#App imports
from shellder import *
from dotjson import dot_to_json

# General app settings
app = Flask(__name__)
SESSION_TYPE = "filesystem"
app.config.from_object(__name__)
Session(app)

FILENAME_REGEXP = r"(?:\.{0,2}\/)?(?:[\w\-]+\/)+(?:[\w\-]+\.\w+)?"

def make_url(string):
    files = re.findall(FILENAME_REGEXP, string)
    res = string
    for f in files:
        method = "get"
        if f.endswith(".dot"):
            method = "render"
        elif f.endswith("/"):
            method = "ls"
        url = "<a href=\"%s?file=%s\">%s</a>" % (method, quote(f), f)
        res = res.replace(f, url)
    #print res
    return res

def format_output(lines):
    res = "<br/>".join(map(make_url, lines))
    #print res
    #res = "<br/>".join(lines)
    return res

env_path = "../../../"
cache_path = "static/cache/"
shellders = dict()

@app.route("/", methods=['GET'])
def index():
    if "username" in session:
        return flask.render_template("index.html", console=format_output(session["log"]))
    else:
        logged = shellders.keys()
        return flask.render_template("login.html", names=logged)

@app.route("/login", methods=['GET'])
def login():
    global shellders
    if "username" not in session:
        name = request.args.get("username", "gaf")
        session["username"] = name
        if name not in shellders:
            shellders[name] = Shellder("/tmp/vis_in_" + name, "/tmp/vis_out_", env_path)
            log = shellders[name].get_output()
            print "Got", log
            session["log"] = log
        else:
            session["log"] = ["(the previous session log has been lost)"]
    return flask.redirect("/")

@app.route("/logout", methods=['GET'])
def logout():
    if "username" in session:
        name = session.pop("username")
        shellder = shellders.pop(name, None)
        shellder.close()
        return "<html><body><p>You have been logged out.</p><p><a href=\"/\">Return to main</a></p></body></html>"
    return "No opened session"

@app.route("/command", methods=['GET'])
def command():
    sh = shellders[session["username"]]
    com = request.args.get("command", "")
    if len(com):
        print "Sending `%s`..." % com
        sh.send(com)
    result = sh.get_output(5)
    complete = False
    if len(result) and result[-1] == sh.end_out:
        complete = True
    session["log"].extend(result)
    #return result
    return flask.jsonify(log=format_output(result), complete=complete)

@app.route("/get")
def get():
    file_path = augment(unquote(request.args.get("file", "")))
    return flask.send_file(file_path, as_attachment=True, attachment_filename=path.basename(file_path))

@app.route("/render")
def render():
    global cache_path
    file_path = unquote(request.args.get("file", ""))
    type = request.args.get("method", "svg")
    _, full_name = path.split(file_path)
    name_only, _ = path.splitext(full_name)
    if type == "svg":
        res_path = cache_path + name_only + ".svg"
        result = open(res_path, "w")
        subprocess.call(["dot", "-Tsvg", env_path + file_path], stdout=result)
        result.close()
        return res_path
    elif type == "json":
        input = open(env_path + file_path, "r")
        content = dot_to_json(input.read())
        if (request.args.get("result", "") == "response"):
           return content
        res_path = cache_path + name_only + ".json"
        result = open(res_path, "w")
        result.write(content)
        result.close()
        return res_path
    else:
        return "Unknown method"

@app.route("/vertex/<vertex_id>")
def vertex(vertex_id):
    shellder = shellders[session["username"]]
    res_path = next((cache_path + f for f in listdir(cache_path) if f.endswith(vertex_id + "_.svg")), None)
    if res_path is None:
        #Render a new file
        shellder.send("draw_vertex " + vertex_id)
        file_path = re.finditer(FILENAME_REGEXP, shellder.get_output()[0]).next().group()
        _, full_name = path.split(file_path)
        name_only, _ = path.splitext(full_name)
        res_path = cache_path + name_only + ".svg"
        result = open(res_path, "w")
        subprocess.call(["dot", "-Tsvg", env_path + file_path], stdout=result)
        result.close()
    return flask.redirect(res_path)

def augment(path):
    global env_path
    if path.startswith("/"):
        return path
    else:
        return env_path + path

#@app.route("/folder/<path:path>")
def folder(path):
    print "Getting contents of", path
    global env_path
    files = [path + f for f in os.listdir(augment(path))]
    print files
    return [f for f in files if isfile(augment(f))]

@app.route("/ls")
def ls():
    dir_path = unquote(request.args.get("file", ""))
    return format_output(folder(dir_path))

if __name__ == "__main__":
    app.debug = True
    app.secret_key = "somekey"
    app.run()
