"""
Phase 2 driver: run src/asr/run_month.py for every month under data/raw,
one subprocess per month (bounded memory), skipping months already done.
Writes outputs/asr/run_all_log.txt.
"""
import glob
import json
import os
import subprocess
import sys
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def main():
    extra = sys.argv[1:]
    months = sorted("/".join(p.replace("\\", "/").split("/")[-2:])
                    for p in glob.glob(os.path.join(ROOT, "data", "raw", "*", "*")) if os.path.isdir(p))
    log = open(os.path.join(ROOT, "outputs", "asr", "run_all_log.txt"), "a")
    for m in months:
        out = os.path.join(ROOT, "outputs", "asr", "months", m.replace("/", "_"))
        cj = os.path.join(out, "calib.json")
        if not extra and os.path.exists(cj) and "n_fitted" in json.load(open(cj)) \
                and json.load(open(cj)).get("n_fitted") == json.load(open(cj)).get("n_dayside"):
            log.write(f"skip {m}\n")
            continue
        t0 = time.time()
        r = subprocess.run([sys.executable, os.path.join(ROOT, "src", "asr", "run_month.py"), "--month", m] + extra,
                           capture_output=True, text=True)
        tail = "\n".join(l for l in r.stdout.splitlines() if l.startswith(("fitted", "20")))
        log.write(f"{m}: rc={r.returncode} {time.time()-t0:.0f}s\n{tail}\n")
        if r.returncode != 0:
            log.write(r.stderr[-3000:] + "\n")
        log.flush()


if __name__ == "__main__":
    main()
