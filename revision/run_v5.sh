#!/usr/bin/env bash
# Revision run (v5 = v3 + CLASS ARF + TEMP-binned detector calibration + 8-17 keV proxy band): fit each month once its archive is fully unpacked.
# 2024_12 first: other months borrow its background-template shapes when they lack umbra exposures.
cd "$(dirname "$0")/.."
X=/c/Users/geete/ch2_data/x
OUT=outputs/asr/months_v5
mkdir -p "$OUT"
order() { echo 2024_12; ls $X/.done_* | sed 's/.*ch2_cla_l1_//' | grep -v 2024_12; }
while :; do
  todo=0
  for m in $(order); do
    [ -f "$OUT/$m/fits.parquet" ] || [ -f "$OUT/$m.failed" ] && continue
    todo=1
    echo "$(date +%H:%M) start $m"
    python src/asr/run_month.py --month "${m/_//}" --gamma_zero --response arf --temp_cal --hi_band 8 17 --out_root "$OUT" > "$OUT/$m.log" 2>&1
    rc=$?; [ $rc != 0 ] && touch "$OUT/$m.failed"
    echo "$(date +%H:%M) $m rc=$rc $(grep '^fitted' "$OUT/$m.log")"
  done
  n=$(ls $X/.done_* | wc -l)
  [ "$todo" = 0 ] && [ "$n" -ge 55 ] && break
  [ "$todo" = 0 ] && sleep 120
done
echo V5_DONE
