#!/usr/bin/env python
import subprocess
import sys


def main(args):
    jobid = args[1]

    output = str(
        subprocess.check_output(
            "sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True
        ).strip()
    )

    running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
        print("success")
    elif any(r in output for r in running_status):
        print("running")
    else:
        print("failed")


if __name__ == "__main__":
    sys.exit(main(sys.argv))
