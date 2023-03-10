try:
    import urllib.request as url_lib #Python 3
except ImportError:
    import urllib as url_lib # Python 2
def read_url(url):
    f = url_lib.urlopen(url)
    txt = f.read()
    return txt   