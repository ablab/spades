"""
*** Original blogpost: http://tomerfiliba.com/blog/Hypertext/

Hypertext: "It's Snakes all the Way Down"
=========================================

Who needs ``haml`` or templating languages? Hypertext is an in-language DSL for producing 
perfectly valid (and safe) XHTML directly from Python code! It's a bit magical, that's true,
but it's thread-safe and it makes rendering HTML pages almost as compact as ``haml``.

Creating HTML elements is easy and straightforward ::

    >>> import hypertext as H
    >>> x = H.h1("hello", class_ = "highlight")
    >>> x
    <Element h1(class='highlight'), 1 subelements>

Invoking ``str`` (or via ``print``) on an HTML element renders it textually ::

    >>> print x
    <h1 class="highlight">hello</h1>

Alternatively, if you're not afraid of polluting your namespace ::

    >>> from hypertext import *
    >>> print h1("hello", class_ = "highlight")
    <h1 class="highlight">hello</h1>

And there's a shortcut for specifying the class (a la ``haml``) ::

    >>> print h1.highlight("hello")
    <h1 class="highlight">hello</h1>

Which can be combined for specifying multiple classes... ::

    >>> print h1.highlight.important("hello")
    <h1 class="highlight important">hello</h1>

Element can be nested ::
    >>> print div.content(h1.highlight("hello"))
    <div class="content"><h1 class="highlight">hello</h1></div>

But also using ``with`` statement

    >>> with div.content as root:                            # doctest:+ELLIPSIS
    ...     h1.highlight("hello")
    ...
    ,,,
    >>> print root
    <div class="content"><h1 class="highlight">hello</h1></div>

.. note::
    Don't mind the ``,,,``, it's there to make ``doctest`` happy

The attributes and text of elements can be specified after their construction :: 

    >>> with div.content as root:
    ...     with h1:
    ...         ATTR(class_ = "highlight")
    ...         TEXT("hello")
    ...
    >>> print root
    <div class="content"><h1 class="highlight">hello</h1></div>

Text is HTML-escaped by default, but you can override it (if you trust its source) ::

    >>> with div.content as root:
    ...     dangerous = "<script>alert('oh no!');</script>"
    ...     TEXT("hello %s" % (dangerous,))
    ...     UNESCAPED("hello %s" % (dangerous,))
    ...
    >>> print root
    <div class="content">
      hello &lt;script&gt;alert(&apos;oh no!&apos;);&lt;/script&gt;
      hello <script>alert('oh no!');</script>
    </div>

Putting it all together::

    >>> with html as doc:                                    # doctest:+ELLIPSIS
    ...     with head:
    ...         meta(http_equiv="Content-Type", content="text/html; charset=utf-8")
    ...         title("welcome to my page")
    ...     with body:
    ...         ATTR(id="body")
    ...         with div.content(id = "floop", data_role = "page").pretty:
    ...             TEXT("Hello, my <name> is")
    ...             strong("Bob")
    ...             TEXT("Also known as", em("Robert"))
    ...             UNESCAPED("and I <b>like</b>")
    ...             with ul:
    ...                 for item in ["cats", "rats", "hats"]:
    ...                     li(item)
    ...
    ,,,
    >>> print doc
    <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html xmlns="http://www.w3.org/1999/xhtml">
      <head>
        <meta content="text/html; charset=utf-8" http-equiv="Content-Type"/>
        <title>welcome to my page</title>
      </head>
      <body id="body">
        <div data-role="page" class="content pretty" id="floop">
          Hello, my &lt;name&gt; is
          <strong>Bob</strong>
          Also known as
          <em>Robert</em>
          and I <b>like</b>
          <ul>
            <li>cats</li>
            <li>rats</li>
            <li>hats</li>
          </ul>
        </div>
      </body>
    </html>

The use of ``with`` and the fact you're dealing with pure Python code makes things easy - 
there's no need for templating languages, passing parameters back and forth, extending base 
templates or including files. It's snakes all the way down!
"""
import threading


__all__ = ["Element", "TEXT", "UNESCAPED", "ATTR", "EMBED", "THIS", "PARENT"]

class Unescaped(str):
    __slots__ = ()
    def __repr__(self):
        return "Unescaped(%s)" % (str.__repr__(self))

_MAPPING = {"&" : "&amp;", "'" : "&apos;", '"' : "&quot;", "<" : "&lt;", ">" : "&gt;"}

def xml_escape(text):
    if isinstance(text, Unescaped):
        return str(text)
    else:
        return "".join(_MAPPING.get(ch, ch) for ch in str(text))

#===================================================================================================
# <magic>
#===================================================================================================
_per_thread = threading.local()

class Element(object):
    class __metaclass__(type):
        __slots__ = ()
        def __getattr__(cls, name):
            return cls(class_ = name)
        def __enter__(cls):
            return cls().__enter__()
        def __exit__(cls, t, v, tb):
            _per_thread._stack[-1].__exit__(t, v, tb)

    __slots__ = ["_attrs", "_elems"]
    TAG = None
    INLINE = False
    DOCTYPE = None
    DEFAULT_ATTRS = {}
    
    def __init__(self, *elems, **attrs):
        self._attrs = self.DEFAULT_ATTRS.copy()
        self._elems = []
        self._parent = None
        if getattr(_per_thread, "_stack", None):
            _per_thread._stack[-1]._elems.append(self)
            self._parent = _per_thread._stack[-1]
        else:
            self._parent = None
        self(*elems, **attrs)
    
    def __repr__(self):
        return "<Element %s(%s), %s subelements>" % (self.TAG if self.TAG else self.__class__.__name__.lower(),
            ", ".join("%s=%r" % (k, v) for k, v in self._attrs.items()), len(self._elems))
    def __str__(self):
        return render(self)
    
    def __enter__(self):
        if not hasattr(_per_thread, "_stack"):
            _per_thread._stack = []
        _per_thread._stack.append(self)
        return self
    def __exit__(self, t, v, tb):
        _per_thread._stack.pop(-1)
    
    def __getitem__(self, name):
        return self._attrs[name]
    def __delitem__(self, name):
        del self._attrs[name]
    def __setitem__(self, name, value):
        self._attrs[name] = value
    
    def __call__(self, *elems, **attrs):
        for elem in elems:
            if isinstance(elem, Element):
                if elem._parent is not None:
                    elem._parent._elems.remove(elem)
                elem._parent = self
            self._elems.append(elem)
        for k, v in attrs.items():
            if k.endswith("_"):
                k = k[:-1]
            self._attrs[k.replace("_", "-")] = v
        return self
    
    def __getattr__(self, name):
        if "class" in self._attrs:
            self._attrs["class"] += " " + name
        else:
            self._attrs["class"] = name
        return self

def THIS():
    """Return the current HTML element on the stack; raises ``IndexError`` in case the stack 
    is empty"""
    return _per_thread._stack[-1]
def PARENT(count = 1):
    """Returns the ``count``-parent of the current HTML element on the stack; raises ``IndexError``
    if there's no ``count`` parent. ``count`` defaults to 1, which means the immedaite parent"""
    return _per_thread._stack[-1 - count]
def TEXT(*texts):
    """Appends the given texts (as well as HTML elements) to the current element"""
    THIS()(*texts)
def UNESCAPED(*texts):
    """Appends the given texts unescaped to the current element. Note: security risk, use 
    with caution"""
    THIS()(*(Unescaped(text) for text in texts))
def EMBED(element):
    """Embeds (appends) the given HTML element into the current element (transfers ownership)"""
    THIS()(element)
def ATTR(**kwargs):
    """Sets the given keyword-arguments as attributes of the current HTML element"""
    THIS()(**kwargs)
#===================================================================================================
# </magic>
#===================================================================================================

def _render(elem, level, dont_indent = False):
    indent = "  " * level
    if not isinstance(elem, Element):
        result = xml_escape(elem)
    else:
        tag = elem.TAG if elem.TAG else elem.__class__.__name__.lower().rstrip("_")
        attrs = " ".join('%s="%s"' % (k, xml_escape(str(v))) for k, v in elem._attrs.items())
        if attrs:
            attrs = " " + attrs
        if elem._elems:
            elements = "\n".join(_render(e, level + 1) for e in elem._elems)
            if elem.INLINE or not "\n" in elements:
                elements = _render(elem._elems[0], level + 1, True)
                result = "<%s%s>%s</%s>" % (tag, attrs, elements, tag)
            else:
                result = "<%s%s>\n%s\n%s</%s>" % (tag, attrs, elements, indent, tag)
        else:
            result = "<%s%s/>" % (tag, attrs)
    if dont_indent:
        return result
    else:
        return indent + result

def render(root):
    """renders the given HTML element and prepends ``DOCTYPE`` if one exists"""
    raw = _render(root, 0)
    if root.DOCTYPE:
        raw = root.DOCTYPE + "\n" + raw
    return raw

#===================================================================================================
# HEAD elements
#===================================================================================================
class html(Element): 
    DEFAULT_ATTRS = {"xmlns" : "http://www.w3.org/1999/xhtml"}
    DOCTYPE = ('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" '
        '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">')

class html5(Element):
    DOCTYPE = '<!DOCTYPE html>'
    TAG = 'html'

class head(Element): pass
class link(Element): pass
class meta(Element): pass
class title(Element): INLINE = True
class script(Element): DEFAULT_ATTRS = {"type" : "text/javascript"}
class style(Element): DEFAULT_ATTRS = {"type" : "text/css"}

#===================================================================================================
# BODY elements
#===================================================================================================
class body(Element): pass
class p(Element): pass
class div(Element): pass
class blockquote(Element): pass
class dl(Element): pass
class dt(Element): pass
class dd(Element): pass
class li(Element): pass
class ul(Element): pass
class ol(Element): pass

class form(Element): pass
class input(Element): pass
class button(Element): pass
class select(Element): pass
class label(Element): pass
class optgroup(Element): pass
class option(Element): pass
class textarea(Element): pass
class legend(Element): pass

class table(Element): pass
class tr(Element): pass
class th(Element): pass
class td(Element): pass
class colgroup(Element): pass
class thead(Element): pass
class tbody(Element): pass
class tfoot(Element): pass

class frame(Element): pass
class iframe(Element): pass
class noframe(Element): pass
class frameset(Element): pass

class pre(Element): INLINE = True
class code(Element): INLINE = True
class span(Element): INLINE = True
class a(Element): INLINE = True
class br(Element): INLINE = True
class hr(Element): INLINE = True
class em(Element): INLINE = True
class strong(Element): INLINE = True
class cite(Element): INLINE = True
class h1(Element): INLINE = True
class h2(Element): INLINE = True
class h3(Element): INLINE = True
class h4(Element): INLINE = True
class h5(Element): INLINE = True
class h6(Element): INLINE = True
class i(Element): INLINE = True
class b(Element): INLINE = True
class u(Element): INLINE = True
class sub(Element): INLINE = True
class sup(Element): INLINE = True
class big(Element): INLINE = True
class small(Element): INLINE = True
class img(Element): INLINE = True

__all__.extend(cls.__name__ for cls in Element.__subclasses__())

if __name__ == "__main__":
    import doctest
    # nice trick from http://stackoverflow.com/a/9400829/434796
    doctest.ELLIPSIS_MARKER = ",,,"
    doctest.testmod()

