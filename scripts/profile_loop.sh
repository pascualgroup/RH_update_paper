#!/usr/bin/env bash

df_prefix="param_grids/Surat/Max_Temp/param_grid_"
if [ -z "$RUN_ID" ]; then
    # generate timestamp + short random hex (openssl if available, else date+pid)
    TS=$(date -u +%Y%m%dT%H%M%SZ)
    if command -v openssl >/dev/null 2>&1; then
        SUF=$(openssl rand -hex 2)
    else
        SUF=$((RANDOM % 10000))
    fi
    RUN_ID="${TS}_${SUF}"
fi
export RUN_ID




for param in sigOBS sigPRO muS2S1 muEI1 muI1S2 muI2S2 betaOUT rho tau q0 alpha b1 b2 b3 b4 b5 b6 bH S1_0 E_0 I1_0 I2_0 S2_0 K_0 F_0; do
    bash scripts/submit_profile.sh $param "${df_prefix}${param}.csv"
done
