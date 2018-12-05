---
author: Benjamin ADOR
title: Preprocessing notebook
---

09/02/2018
==========
### Epoch rejection: mag = 2500 fT; ICA rejection: mag = 4000 fT; ECG threshold = 0.3; EOG threshold = 3.0

ECG components: epochs centered around R (`qrs_detector`), P around -200 ms, T around +300 ms

Low pulse (091: 50 bpm <=> 1.2 s heartbeat window)
  * R events: 
  * T events:
        -400 -- +800 ms from R peak to P wave
        -550 -- +600 ms from P wave to before the next one

High pulse (050/069: 96 bpm <=> 0.625 s heartbeat window)
  * R events: 
  * T events:
        -250 -- +400 ms from R peak to P wave
        -400 -- +300 ms from P wave to before the next one

**T epoch window**:
    -550 ms takes the whole cardiac artefact
    -175 ms never takes R peak
    +275 ms never takes the following cardiac artefact

**R epoch window**:
    -350 ms never takes the previous cardiac artefact
    +425 ms never takes the following cardiac artefact


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


Subject 011
-----------
_Rejected before I arrived, don't know why_


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
`OM_03`: only 3 EOG events --> actually None (or one the 1st second). EOG100 (2nd) is not EOG

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
`FA_02`: only 8 EOG events --> actually 4 (3 ECG components, ECG36 seems ok)


### EOG bug => re-run from subject 037

037
---
### ECG artifacts not removed (see overlay)


042
---
(`FA_06`: 3 ECG components, ECG1 seems ok)


050
---
Pulse > 90bpm (checked, OK)

`FA_02`: not a single blink


052
---
Nothing to report


053
---
>100 epochs dropped: 0.5 Hz highpass not enough

>3 ECG components: ECG waves

`FA_02`: only 7 blinks (checked)

`OM_06`: only 3 blinks (checked)


054
---
>3 ECG components: ECG waves

`FA_11`: weak EOG18 -> 7 actual blinks, not very sharp


055
---
Nothing to report


063
---
`OM_06`: ECG rejection not very effective

`OM_14`: no ECG component -> one at 0.175 threshold (no effect), three at 0.15 (small effect, probably not worth it)


064
---
Nothing to report


068
---
Nothing to report


069
---
Pulse > 90bpm (checked, OK)

ECG waves (theta?)

weak EOG components --> wavy channel

high amplitude alpha

`RS_01`: 3 ECG components


070
---
2 EOG, 3 ECG components: OK


073
---
very few EOG events (`FA_06`: only 2 actual blinks)

`OM_02`: ECG10 (3rd components) weak and weird

`FA_11`: >3 ECG events, ECG47 weak and weird


081
---
**Rejected** (temporal artifact)


083
---
>3 ECG components (seems OK)

`FA_06`: 2nd EOG component (EOG2) for 3 events only


087
---
Nothing to report


089
---
`RS_02`: >3 ECG components, seems OK

`OM_11` : 3 ECG components, ECG36 is weak and weird


090
---
>3 ECG components --> ECG waves, noisy components

`FA_10`: >100 epochs dropped --> bad channel 'MLT51-2805' (would also save 37 out of 87 'bad' epochs for `RS_01`)


091
---
Nothing to report


093
---
>3 ECG components --> ECG waves


094
---
weak ECG components, ineffective rejection

>100 epochs dropped --> right temporal muscle (=> lowpass?)


096
---
weak ECG components corresponding to ECG noise

`RS_02`: >100 epochs dropped --> high-amplitude noise


23/02/2018
==========
Finally done checking \o/

Overall, ineffective ICA comes from weak components, decreasing thresholds would not improve artifact detection.

`014/OM_13-EOG100` should not be detected

'ECG waves': 053, 054, 069, 071, 090, 093

Bad data: 037, 053, 069, 090, 094, 096


26/02/2018
==========

056
---
3rd ECG component: some ECG waves


067
---
Nothing to report


071
---
>3 ECG components --> ECG waves


072
---
`FA_16`: ECG rejection not very effective (T component in particular); weird "split" EOG component


095
---
No ECG --> use EEG008-2800 instead (although EEG003 has more ECG amplitude)

`FA_02`: >3 ECG components --> maybe some weird T waves


109
---
`RS_01`, `FA_02`: >3 ECG components (OK)

`OM_14`: 2nd weird EOG component


01/03/2018
==========

004
---
`RS_01`: heartbeat count by eye (420 s): 545 --> 78 bpm

Flat ECG --> epoching either on S peak (2/3) or on T peak (1/3) instead of R.


010
---
`FA_03`: 2 almost simultaneous ECG R events (=bug), leading to twice the same T peak --> not unique events error

This is because the ECG S peak is bigger than the R peak --> wrong epoching


07/03/2018
==========

054
---
`OM_15`: right ear moves >0.5cm from the beginning to the end of the block.


15/03/2018
==========

010
---
01: Huge movement at the last minute

07: Couldn't keep coil proportions (MinMax, difference <= 0.0542 cm)


018
---
07: Huge movement during the last minute


028
---
14: Couldn't keep coil proportions (MinMax, difference <= 0.0519 cm); still one bad segment at 5min38


030
---
01: Couldn't keep coil proportions (MinMax, difference <= 0.0128 cm)


040
---
02: Couldn't keep coil proportions (sharp peak => Mean, difference<= 0.0175)


042
---
14: Couldn't keep coil proportions (MinMax on Z + MidAll on X and Y, difference <= 0.0276 cm (close to Mean: 0.0296))


### I just realised coil proportions don't work as intended because of the difference between DATA and HC ref...

_**Corrected when coil proportions shouldn't be kept; all others should be zero (but are not corrected)**_


054
---
15: Too much drift --> Couldn't keep coil proportions nor prevent bad segments

`+end`: preserved the end to the cost of the first 20 sec (Mean, difference <= 0.0135 cm)

`+start`: preserved the beginning (except the first 1.2 sec) to the cost of the last 14 sec (MinMax, difference <= 0.0095 cm)


068
---
14: Tight fit --> MidAll on Z + 0.15 on X + 0.025 on Y


072
---
01: Couldn't keep coil proportions (Median, difference <= 0.0561 cm)

`03+end`: preserved the end to the cost of the first 100 sec (Median, difference <= 0.0514 cm)

`03+start`: preserved the beginning to the cost of the last 42 sec (MinMax with Z_Left set to -0.45, difference <= 0.0301 cm)


073
---
01: MinMax (difference <= 0.0292), bad = [60,90] + last 2 min

`01+end`: preserved the end to the cost of the first 120 sec (MinMax with X_Left = 0.25, Y_Left = 0.2, Z_Left = -0.4, difference <= 0.0477 cm)

`01+start`: preserved the beginning to the cost of the last 140 sec (MinMax with X_Left = 0.15, Y_Left = 0.1, Z_Left = -0.2, difference <= 0.0404 cm)


21/03/2018
==========

032
---
`RS_01`: not best LPA/RPA fit


067
---
Surface inner skull is not completely inside surface outer skull => fixed with make_bem_model(..., conductivity=[0.3])


087
---
Weird head


073
---
No head-dense.fif, coregistration gui works anyway (same for subjects 083 and 095)


23/03/2018
==========

054
---
Noise data are from the day before (same for subjects 068 and 109) => updated csv_io not to get the error


02/05/2018
==========

057
---
`14`: couldn't keep coil proportions (MinMax, difference <= 0.0177 cm)


057
---
`01`: couldn't keep coil proportions (MinMax, difference <= 0.0303 cm)


065
---
`06`: Huge movement at the end of meditation, impossible to preserve the end (about 10 sec lost, Mean, difference <= 0.0320 cm)


074
---
Weird oscillatory movement (particularly Nasion)


075
---
`01`: couldn't keep coil proportions (MinMax, difference <= 0.0371 cm)

`06`: couldn't keep coil proportions (MinMax, difference <= 0.0205 cm)


098
---
01: Huge movement at 75 sec, and another substantial one at 285 sec.

`01+end`: preserved the end to the cost of the first 2 min (Median with Z_Left = -0.55, difference <= 0.0535 cm)

`01+start`: preserved the beginning to the cost of the last 2 min + 5 sec around 75 sec (Median, difference <= 0.0235 cm)

`06+end`: preserved the end to the cost of the first 45 sec (Mean, difference <= 0.0492 cm)

`06+start`: preserved the beginning to the cost of the last 80 sec (MidAll with Z = -0.44)


103
---
`06+start`: Peak movement at 363 sec (one bad segment)

`06+end`: Preserved the bad segment but couldn't keep coil proportions (MinMax with X_Nasion = 0 and Z_Nazion = 0.32, difference <= 0.0952 cm)


105
---
`01`: couldn't keep coil proportions (MinMax, difference <= 0.0121 cm)


04/05/2018
==========

057
---
No blinks ? Or very small amplitude --> EOG component picks up alpha


058
---
`FA_14`: 2nd weak and noisy EOG component


059
---
Nothing to report.


065
---
Weak ECG components (None for `FA_10`).

`OM_14`: Weak EOG component (no obvious blink in data, but still ocular activity).


092
---
Nothing to report.


101
---
Nothing to report.


104
---
Nothing to report.


108
---
`FA_14`: weird EOG components (probably due to numerous double blinks + long blinks)


02/10/2018
==========

079
---
`01`: 3 peak movements around 52, 54, and 65 sec ; couldn't keep coil proportions (Mean, difference <= 0.0346 cm)


080
---
`10`: couldn't keep coil proportions (Median, difference <= 0.0189 cm)


099
---
`10+start`: preserved the beginning to the cost of the last 100 sec (X = 0.2, Y = 0.05, Z = -0.4)

`10+end`: Preserved the end to the cost of the first 145 sec (MinMax with X_Left = 0.45, Y_Left = 0.1 and Z_Left = -0.9, difference <= 0.0214 cm)


21/11/2018
==========
ECG Reports

004
---
The QRS complex is weak and wavy, such that peak detection is not consistent even given the R sign a priori.

=> rejected


010
---
T peak seems to have a polarity opposite to that of the P and R peaks.

Wrong peak detection (apparently from FA04), try brute force.

=> rejected


012
---
Inconsistent peak detection => Actually no, individual ECG events displays this weird fork. It looks like, instead of having a PQRST polarity of +-+-+ (or -+-+-), the ECG does -++-+, as if the polarity changes between Q and R.


028
---
Detection a bit off-peak, but consistent.


056
---
Visually spotted 12 inconsistent R peaks for `RS01` (at times 5, 86, 202, 213, 333, 339.5, 357, 365, 371, 396, 400, 417), all corresponding to rejected T events (thus equalize_event_counts would do the trick).

The artefact this inconsistency causes on the cardiac ERP is also present for the other blocks.

=> specifying the R sign should fix this. (Fixed)


064
---
The S peak has a greater amplitude than the R peak and is thus detected instead. => Re-run with custom_args. (Fixed)


069
---
One spurious R event at the start of `FA02`. (Re-run)


076
---
The S peak has a greater amplitude than the R peak and is thus detected instead. => Re-run with custom_args. (Fixed)


094
---
High amplitude S peak, specify R sign. (Fixed)
