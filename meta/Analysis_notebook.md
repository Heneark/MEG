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
**Pulse** 80-90 bpm (FA_04: at sight it seems to be 78, while ICA says 91)

`FA_02`: 2nd weak and noisy **EOG** component

`OM_04`: 2nd weak and noisy **EOG** component


Subject 010
-----------
**drop_log**: ~200 epochs dropped

`RS_01`: 2nd weak **EOG** component, around +200 ms

`OM_07`: **ECG**25 has a high score but not marked for exclusion => T component


Subject 12
----------
