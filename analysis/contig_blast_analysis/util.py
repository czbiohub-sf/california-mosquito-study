#!/usr/bin/env python3
# script from https://github.com/czbiohub/pyblastc/blob/master/scripts/util.py

import sys
import json
import multiprocessing
import time

def parse_headerless_table(rows, schema={}):
    headers = schema.keys()
    for values in rows:
        assert len(headers) == len(values)
        yield {h: schema[h](v) for h, v in zip(headers, values)}

def tsv_rows(path):
    with open(path, "r") as stream:
        for line in stream:
            if not line.startswith("#"):
                yield line.rstrip("\n").split("\t")

# Thread-safe and timestamped prints.
tslock = multiprocessing.RLock()


def timestamp(t):
    # We do not use "{:.3f}".format(time.time()) because its result may be
    # up to 0.5 ms in the future due to rounding.  We want flooring here.
    s = str(int(t * 10))
    return s[:-1] + "." + s[-1:]


def tsfmt(msg):
    ts = timestamp(time.time()) + " "
    msg = ts + msg.replace("\n", "\n" + ts)
    return msg


def tsout(msg):
    with tslock:
        sys.stdout.write(str(msg))
        sys.stdout.write("\n")


def tserr(msg):
    with tslock:
        sys.stderr.write(str(msg))
        sys.stderr.write("\n")


def tsprint(msg):
    tserr(tsfmt(msg))

if __name__ == "__main__":
    tsprint("Hello from pyrecomb utils.")

def check_exit_code(process, command):
    """ Capture stdout, stderr. Check unix exit code and exit if non-zero """
    out, err = process.communicate()
    if process.returncode != 0:
        err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
        sys.exit(err_message)
