import urllib3
import io

def write_html(url, path):
    http = urllib3.PoolManager()
    r = http.request("Get", url=url, preload_content=False)
    r.auto_close = False
    with open(path, 'w') as f:
        for line in io.TextIOWrapper(r):
            f.write(line)


def test():
    url = "http://localhost:8888/tree.php"
    path = './test.html'
    write_html(url, path)


if __name__ == '__main__':
    test()