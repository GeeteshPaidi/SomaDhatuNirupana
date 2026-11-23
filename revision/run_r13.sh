#!/usr/bin/env bash
# R1-3: refit three months with exposure-matched 1 s GOES temperatures, after the v4 run.
cd "$(dirname "$0")/.."
until grep -q V5_DONE outputs/asr/run_v5.log; do sleep 120; done
for m in 2024_12 2025_06 2025_12; do
  [ -f outputs/asr/solar/goes_tem_1s_$m.parquet ] || \
    python src/asr/goes_1s.py --month $m --fits outputs/asr/months_v5/$m/fits.parquet 2>&1 | grep -v Warn
  python src/asr/run_month.py --month ${m/_//} --gamma_zero --response arf --temp_cal --hi_band 8 17 \
    --goes_tem outputs/asr/solar/goes_tem_1s_$m.parquet --out_root outputs/asr/months_r13 \
    > outputs/asr/months_r13_$m.log 2>&1
  echo "$(date +%H:%M) $m rc=$? $(grep '^fitted' outputs/asr/months_r13_$m.log)"
done
echo R13_DONE
