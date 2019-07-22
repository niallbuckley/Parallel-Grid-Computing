#!/usr/bin/env python3
import threading, time, socket
from sys import argv, stdout
from itertools import chain
from concurrent.futures import ThreadPoolExecutor, as_completed

def main():
    if len(argv) > 1:
        out_stream = open(argv[1], "w")
    else:
        out_stream = stdout
    machines = find_machines()
    out_stream.writelines('\n'.join(machines) + '\n')
    if out_stream != stdout:
        out_stream.close()


def find_machines():
    csg24 = ('csg24-{:02d}'.format(n) for n in range(2, 48))
    cs109 = ('cs1-09-{:02d}'.format(n) for n in range(1, 25))
    cs131 = ('cs1-31-{:02d}'.format(n) for n in range(14, 17))
    hostnames  = ['{}.ucc.ie'.format(p) for p in chain(csg24, cs109, cs131)]
    unavailable, available = check_hosts(hostnames)
    return sorted(available)

def check_hosts(hostnames):
    available_hostnames = []
    unavailable_hostnames = []

    with ThreadPoolExecutor(max_workers=len(hostnames)) as executor:
        future_to_name = {executor.submit(check_host, name): name for name in hostnames}
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            if future.result():
                available_hostnames.append(name)
            else:
                unavailable_hostnames.append(name)
    return (unavailable_hostnames, available_hostnames)

def check_host(host_name):
    s = socket.socket()
    try:
        s.connect((host_name, 22))
        s.close()
        return True
    except Exception as e:
        return False

if __name__ == '__main__':
    main()