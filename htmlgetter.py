import argparse
import urllib3
import io

def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--url', '-u', type=str, required=True, dest='url',
        help='Required; Specify a link.')
    parser.add_argument(
        '--out', '-o', type=str, required=True, dest='out',
        help='Required; Specify a output path.')
    args = parser.parse_args()
    return args


def write_html(url, path):
    http = urllib3.PoolManager()
    r = http.request("Get", url=url, preload_content=False)
    r.auto_close = False
    with open(path, 'w') as f:
        for line in io.TextIOWrapper(r):
            f.write(line)

def main():
    args = argument_parser()
    url = args.url
    path = args.out
    write_html(url, path)


def test():
    url = "http://localhost:8888/tree.php"
    path = './php_server/tree.html'
    write_html(url, path)


if __name__ == '__main__':
    main()
    # test()
