#!/usr/bin/env bash
set -e
qcc -Wall -O2 -disable-dimensions BpostBview.c -o BpostBview -L"$BASILISK/gl" -lglutils -lfb_tiny -lm
./BpostBview tb0.00 ts0.01 te1.04


