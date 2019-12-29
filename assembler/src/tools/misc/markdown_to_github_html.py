import requests
import sys
import codecs 
from traceback import print_exc

# GitHub style CSS for code sections
HEADER = '<html> \
<head> \
    <title>SPAdes 3.13.1 Manual</title> \
    <style type="text/css"> \
        .highlight  pre { \
          background-color: #f0f2f4; \
          border-radius: 2px; \
          font-size: 100%; \
          line-height: 1.45; \
          overflow: auto; \
          padding: 16px; \
        } \
    </style> \
</head> \
<body>'

# Extra spaces in the end
FOOTER = '<br/><br/><br/><br/><br/> \
</body> \
</html>'

# Make a POST request to GitHub API
def make_github_request(in_file_name):
    url = 'https://api.github.com/markdown'
    inf = open(in_file_name)
    mdtext = inf.read()
    inf.close()

    convert_obj = {
      "text": mdtext,
      "mode": "gfm",
      "context": "github/gollum"
    }

    x = requests.post(url, json = convert_obj, headers = {"Accept" : "application/vnd.github.VERSION.full+json"})
    return x.text

# Substitute conversion artifacts and dump to file (non-ascii)
def write_html_to_file(txt, out_file_name):
    final_html = HEADER + txt.replace("<br><br>", "<br>").replace('<a name="user-content-', '<a name="').replace('<a href="assembler/', '<a href="').replace('<p><strong>SPAdes', ' <meta charset="UTF-8"> <p><strong>SPAdes') + FOOTER
    outf = codecs.open(sys.argv[2], 'w', encoding='utf-8')
    outf.write(final_html)
    outf.close()


def main():
    if len(sys.argv) != 3:
        print("Usage: " + sys.argv[0] + " <input.md> <output.html>")
        exit(0)

    write_html_to_file(make_github_request(sys.argv[1]), sys.argv[2])


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
