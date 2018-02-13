09/02/2018

Preprocessing
=============

ECG components: epochs centered around R (`qrs_detector`), P around -200 ms, T around +300 ms


Subject 002
-----------
Pulse: 80 bpm (OK, checked)

`RS_01`: 2nd weak **EOG** component, around +200 ms

`FA_02`: same


Subject 004
-----------
**Pulse** 80-90 bpm (FA_04: at sight it seems to be 78, while ICA says 91)

`RS_01`: no **ECG** component marked for exclusion (23 + 34 have a high score though)

`OM_02`: same with 30 + 32

`FA_04`: same with 23 + 36


Subject 007
-----------
`FA_02`: 2nd weak and noisy **EOG** component

`OM_04`: 2nd weak and noisy **EOG** component


Subject 010
-----------
**drop log**: ~200 epochs dropped => temporal high-frequency noise (`OM_07`: MRT32,41,21,31,22,42; MLT32,42 -2805)

`RS_01`: 2nd weak **EOG** component, around +200 ms

`OM_07`: **ECG**25 has a high score but not marked for exclusion => T component


Subject 012
-----------
`RS_01`: 2nd weak and noisy **EOG** component

`OM_03`: _impressive_ overlay...; 2nd weak and noisy **EOG** component

`OM_04`: 2nd weak and noisy **EOG** component

`FA_06`: noisy overlay; 2nd weak and noisy **EOG** component

`FA_07`: 2nd weak and noisy **EOG** component


Subject 014
-----------
`OM_03`: only 17 blinks => 2 additional weak **EOG** components

`OM_O4`: same with 43 blinks

`FA_O6`: 68 blinks, 1 additional weak **EOG** component


Subject 016
-----------
`RS_01`: some slow drift (LT, MRO53, MRT57); 2nd weak and noisy **EOG** component

`OM_02`: unstable trend (MLT in particular); 2nd weird **EOG** component (temporal?)

`FA_04`: 2nd **EOG** component based on a single temporal event


Subject 018
-----------
`OM_03`: noise persisting on overlay

`OM_04`: same

`FA_06`: weird **EOG** events

`FA_07`: 4 **EOG** components


Subject 028
-----------
`RS_01`: 3 **EOG** components, shifted in time

`FA_02`: strange ECG component, 4 **EOG** components

`OM_06`: 144 **epochs** dropped (NOISE, almost 1 min unexploitable + some moderately bad channels); 2nd weak and noisy **EOG** component

`FA_10`: 145 **epochs** dropped (same)

`OM_14`: 183 **epochs** dropped (same); 2nd weak and noisy **EOG** component


Subject 030
-----------
`FA_02`: 3rd weird **ECG** component; 2nd weird **EOG** component, not frontal, +150 ms

`OM_06`: 130 **epochs** dropped (lot of noise); 2nd weird **ECG** component

`OM_10`: noisy overlay (71 epochs dropped); 2nd weird **ECG** component

`FA_14`: two additional weak **ECG** components


Subject 032
-----------
**Nothing to report**


Subject 037
-----------
**drop log** 200-250 epochs dropped


Subject 40
----------
**NOISE** => rejected

`RS_01`: MLF11->14; MLT11,21,31,32,41,42,43,51->57
