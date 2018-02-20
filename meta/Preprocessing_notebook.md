---
author: Benjamin ADOR
title: Preprocessing notebook
---

09/02/2018
==========
### Epoch rejection: mag = 2500 fT; ICA rejection: mag = 4000 fT; ECG threshold = 0.3; EOG threshold = 3.0

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


14/02/2018
==========
### Epoch rejection: mag = 3500 fT; ICA rejection: mag = 7000 fT; ECG threshold = 0.25; EOG threshold = 3.5

Subjects 004 (all), 069 (`RS_01` and `FA_02`) and 094 (`OM_02`) have no ECG components detected.

Subjects 053 (`FA_02`and `OM_06`), 069 (`FA_02`) and 073 (`RS_01`) have no EOG components detected.

Possible solutions?

1. Low threshold but only take n components max (e.g. 2 for ECG)
2. High threshold but select the highest component if none has been selected (should work for EOG)


Subject 004
-----------
Same issue as with 0.3 ECG threshold: 2 ECG components with higher score but not detected

=> down to 0.1 for detection
### These 2 components are not ECG (see manually generated component plots).


Subject 018
-----------
High number of EOG components

**TO CHECK**


Subject 037
-----------
With ECG threshold of 0.3, `FA_10`and `OM_14` didn't reach it; now fixed.

`RS_01`:Additional candidate: ECG16 (score 0.21)


Subject 053
-----------
With default EOG threshold of 3.0, `OM_06` didn't reach it.

`FA_02`: 

`OM_06`: 

**TO-DO**: check the highest component


Subject 063
-----------
With ECG threshold of 0.3, `OM_06`and `OM_14` didn't reach it; now fixed.

`OM_06`: Additional candidate: ECG22 (score 0.24)


Subject 069
-----------
Manuallly generated ECG score plots for `RS_01`, `FA_02`, `OM_06` and `FA_14`.

`FA_02`: EOG as well


Subject 070
-----------
High number of EOG components

**TO CHECK**


Subject 073
-----------
`RS_01`: 

**TO-DO**: check the highest component


Subject 094
-----------
`OM_02`: ECG19 has a score of 0.24. Two additional candidates: ECG18 and 21 (score 0.18).

`FA_11`: Additional candidate: ECG20 (score 0.22)


16/02/2018
==========
### Reboot of the storage space => lost previous run

It is unclear how find_bads_eog is scoring / classifying components => EOG_threshold set to 5 for next run, with number of EOG components limited to 1 or 2.


19/02/2018
==========
### Epoch rejection: mag = 3500 fT; ICA rejection: mag = 7000 fT; ECG threshold = 0.2 with max 3 components; EOG threshold = 5 with min 1 component and max 2 components

002
---
Weak 3rd ECG component (`RS_01`: ECG7, `FA_02`: ECG8, `OM_04`: ECG22)


004
---
(Still no ECG)


007
---
ECG artifact is not entirely erased (see overlay).


010
---
`RS_01`, `FA_04`: ECG component T kept


012
---
`RS_01`: 2nd weird EOG component

`OM_03`: weird EOG


014
---
`OM_03`: **only 3 EOG events**

`FA_06`: 2nd EOG component (EOG1) for 2 events only


016
---
(3 ECG components, `OM_02`-ECG25 seems ok)


018
---
Nothing to report


028
---
Weirdly "forked" QRS ECG components

(3 ECG components, `RS_01`-ECG49 seems ok)


030
---
More than 3 ECG components, the weak ones kept are weird.

`RS_01`: 5 ECG components, 2 dropped, ECG3 is weird

`FA_02`: same with ECG1

`OM_06`: 4 ECG components, 1 dropped, ECG0&6 are weird

`OM_10`: same with ECG1&2

`FA_14`: same with ECG2&20


032
---
`FA_02`: **only 8 EOG events** (3 ECG components, ECG36 seems ok)


### EOG bug => re-run from subject 037

037
---
### ECG artifacts not removed (see overlay)


042
---
(`FA_06`: 3 ECG components, ECG1 seems ok)


050
---
**Pulse > 90bpm**


053
---
**>100 epochs dropped**


063
---
`OM_14`: **no ECG component**


069
---
**Pulse > 90bpm**


089
---
**Pulse > 90bpm**


090
---
`FA_10`: **>100 epochs dropped**


094
---
**>100 epochs dropped**


096
---
`RS_02`: **>100 epochs dropped**