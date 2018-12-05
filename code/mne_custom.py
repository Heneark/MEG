#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 14:48:28 2018

@author: benjamin.ador
"""

import mne
from mne.epochs import BaseEpochs, Epochs
from mne.externals.six import string_types
from mne.filter import filter_data
from mne.io import RawArray
from mne.io.base import BaseRaw
from mne.preprocessing.bads import find_outliers
from mne.preprocessing.ctps_ import ctps
from mne.preprocessing.ecg import _get_ecg_channel_index, _make_ecg
from mne.utils import logger, sum_squared, verbose, warn
import numpy as np
from scipy.signal import detrend, hilbert

def qrs_custom(sfreq, ecg, var=1/3, thresh_value='auto', levels=2.5, n_thresh=3,
                 l_freq=8, h_freq=16, filter_length='10s', R_sign=0, T_sign=0, heart_rate=0, tstart=0, force=False):
    """
    Copy-paste from mne.preprocessing.ecg, with a few modifications:
        R_sign: instead of taking abs() of the ecg, provide the expected sign of the R peak to turn it upwards.
        ideal_rate: instead of taking a median heart rate, provide the expected pulse.
        var: provided ideal_rate, reduce the size of the search window and increase the jumps between peaks according to this variability parameter
            (var=2/3 with 80 bpm pulse is equivalent to MNE's default).
    
    
    Detect QRS component in ECG channels.

    QRS is the main wave on the heart beat.

    Parameters
    ----------
    sfreq : float
        Sampling rate
    ecg : array
        ECG signal
    thresh_value : float | str
        qrs detection threshold. Can also be "auto" for automatic
        selection of threshold.
    levels : float
        number of std from mean to include for detection
    n_thresh : int
        max number of crossings
    l_freq : float
        Low pass frequency
    h_freq : float
        High pass frequency
    tstart : float
        Start detection after tstart seconds.
    filter_length : str | int | None
        Number of taps to use for filtering.

    Returns
    -------
    events : array
        Indices of ECG peaks
    """
    if heart_rate:
        period = sfreq * 60 / heart_rate
        win_size = int(round(period * var))
        print('Window size: {}\nJump size: {}'.format(win_size/sfreq, int(round(period - win_size / 2))/sfreq))
    else:
        win_size = int(round((60.0 * sfreq) / 120.0))
    
    if force:
        print('BRUTE FORCE DETECTION BASED ON PRIOR HEART RATE KNOWLEDGE. Please check epoching.')
        data = abs(hilbert(ecg))
        data = filter_data(ecg, sfreq, 5, None, None, filter_length,
                           0.5, 0.5, phase='zero-double', fir_window='hann',
                           fir_design='firwin2')
        time = []
        ii = int(round(tstart *sfreq))
        while ii <= data.size - win_size:
            window = data[ii:ii + win_size]
            max_time = np.argmax(window)
            time.append(ii + max_time)
            ii += int(round(max_time + period - win_size / 2))
        return time
    
    filtecg = filter_data(ecg, sfreq, l_freq, h_freq, None, filter_length,
                          0.5, 0.5, phase='zero-double', fir_window='hann',
                          fir_design='firwin2')

    if R_sign:
        ecg_abs = R_sign * filtecg
    else:
        ecg_abs = np.abs(filtecg)
    init = int(sfreq)

    n_samples_start = int(sfreq * tstart)
    ecg_abs = ecg_abs[n_samples_start:]

    n_points = len(ecg_abs)

    maxpt = np.empty(3)
    maxpt[0] = np.max(ecg_abs[:init])
    maxpt[1] = np.max(ecg_abs[init:init * 2])
    maxpt[2] = np.max(ecg_abs[init * 2:init * 3])

    init_max = np.mean(maxpt)

    if thresh_value == 'auto':
        thresh_runs = np.arange(0.3, 1.1, 0.05)
    elif isinstance(thresh_value, string_types):
        raise ValueError('threshold value must be "auto" or a float')
    else:
        thresh_runs = [thresh_value]

    # Try a few thresholds (or just one)
    clean_events = list()
    for thresh_value in thresh_runs:
        thresh1 = init_max * thresh_value
        numcross = list()
        time = list()
        rms = list()
        ii = 0
        while ii < (n_points - win_size):
            window = ecg_abs[ii:ii + win_size]
            if window[0] > thresh1:
                max_time = np.argmax(window)
                time.append(ii + max_time)
                nx = np.sum(np.diff(((window > thresh1).astype(np.int) ==
                                     1).astype(int)))
                numcross.append(nx)
                rms.append(np.sqrt(sum_squared(window) / window.size))
                ii += int(round(max_time + period - win_size / 2)) if heart_rate else win_size
            else:
                ii += 1

        if len(rms) == 0:
            rms.append(0.0)
            time.append(0.0)
        time = np.array(time)
        rms_mean = np.mean(rms)
        rms_std = np.std(rms)
        rms_thresh = rms_mean + (rms_std * levels)
        b = np.where(rms < rms_thresh)[0]
        a = np.array(numcross)[b]
        ce = time[b[a < n_thresh]]

        ce += n_samples_start
        clean_events.append(ce)

    # pick the best threshold; first get effective heart rates
    rates = np.array([60. * len(cev) / (len(ecg) / float(sfreq))
                      for cev in clean_events])

    # now find heart rates that seem reasonable (infant through adult athlete)
    idx = np.where(np.logical_and(rates <= 160., rates >= 40.))[0]
    if not heart_rate:
        if len(idx) > 0:
            heart_rate = np.median(rates[idx])  # get close to the median
        else:
            heart_rate = 80.  # get close to a reasonable default
    idx = np.argmin(np.abs(rates - heart_rate))
    clean_events = clean_events[idx]
    return clean_events


def custom_ecg_events(raw, custom_args, event_id=999, ch_name=None,
                    l_freq=8, h_freq=16, qrs_threshold='auto',
                    filter_length='10s', return_ecg=False, verbose=None):
    """
    Copy-paste from mne.preprocessing.ecg, with a few modifications:
        Use qrs_custom() instead of qrs_detector()
    
    
    Find ECG peaks.

    Parameters
    ----------
    raw : instance of Raw
        The raw data
    event_id : int
        The index to assign to found events
    ch_name : None | str
        The name of the channel to use for ECG peak detection.
        If None (default), a synthetic ECG channel is created from
        cross channel average. Synthetic channel can only be created from
        'meg' channels.
    tstart : float
        Start detection after tstart seconds. Useful when beginning
        of run is noisy.
    l_freq : float
        Low pass frequency.
    h_freq : float
        High pass frequency.
    qrs_threshold : float | str
        Between 0 and 1. qrs detection threshold. Can also be "auto" to
        automatically choose the threshold that generates a reasonable
        number of heartbeats (40-160 beats / min).
    filter_length : str | int | None
        Number of taps to use for filtering.
    return_ecg : bool
        Return ecg channel if synthesized. Defaults to False. If True and
        and ecg exists this will yield None.
    verbose : bool, str, int, or None
        If not None, override default verbose level (see :func:`mne.verbose`
        and :ref:`Logging documentation <tut_logging>` for more).

    Returns
    -------
    ecg_events : array
        Events.
    ch_ecg : string
        Name of channel used.
    average_pulse : float
        Estimated average pulse.
    """
    idx_ecg = _get_ecg_channel_index(ch_name, raw)
    if idx_ecg is not None:
        logger.info('Using channel %s to identify heart beats.'
                    % raw.ch_names[idx_ecg])
        ecg, times = raw[idx_ecg, :]
    else:
        ecg, times = _make_ecg(raw, None, None, verbose=verbose)

    # detecting QRS and generating event file
    ecg_events = qrs_custom(raw.info['sfreq'], ecg.ravel(),
                            thresh_value=qrs_threshold, l_freq=l_freq,
                            h_freq=h_freq, filter_length=filter_length, **custom_args)

    n_events = len(ecg_events)
    average_pulse = n_events * 60.0 / (times[-1] - times[0])
    logger.info("Number of ECG events detected : %d (average pulse %d / "
                "min.)" % (n_events, average_pulse))

    ecg_events = np.array([ecg_events + raw.first_samp,
                           np.zeros(n_events, int),
                           event_id * np.ones(n_events, int)]).T
    out = (ecg_events, idx_ecg, average_pulse)
    if return_ecg:
        out += (ecg,)
    return out


def custom_ecg_epochs(raw, custom_args, ch_name=None, event_id=999, picks=None, tmin=-0.5,
                      tmax=0.5, l_freq=8, h_freq=16, reject=None, flat=None,
                      baseline=None, preload=True, keep_ecg=False,
                      reject_by_annotation=True, verbose=None):
    """
    Copy-paste from mne.preprocessing.ecg, with a few modifications:
        Uses qrs_custom() instead of qrs_detector()
        
        Returns pulse
    
    Conveniently generate epochs around ECG artifact events.

    Parameters
    ----------
    raw : instance of Raw
        The raw data
    ch_name : None | str
        The name of the channel to use for ECG peak detection.
        If None (default), ECG channel is used if present. If None and no
        ECG channel is present, a synthetic ECG channel is created from
        cross channel average. Synthetic channel can only be created from
        'meg' channels.
    event_id : int
        The index to assign to found events
    picks : array-like of int | None (default)
        Indices of channels to include. If None, all channels are used.
    tmin : float
        Start time before event.
    tmax : float
        End time after event.
    l_freq : float
        Low pass frequency.
    h_freq : float
        High pass frequency.
    reject : dict | None
        Rejection parameters based on peak-to-peak amplitude.
        Valid keys are 'grad' | 'mag' | 'eeg' | 'eog' | 'ecg'.
        If reject is None then no rejection is done. Example::

            reject = dict(grad=4000e-13, # T / m (gradiometers)
                          mag=4e-12, # T (magnetometers)
                          eeg=40e-6, # V (EEG channels)
                          eog=250e-6 # V (EOG channels)
                          )

    flat : dict | None
        Rejection parameters based on flatness of signal.
        Valid keys are 'grad' | 'mag' | 'eeg' | 'eog' | 'ecg', and values
        are floats that set the minimum acceptable peak-to-peak amplitude.
        If flat is None then no rejection is done.
    baseline : tuple | list of length 2 | None
        The time interval to apply rescaling / baseline correction.
        If None do not apply it. If baseline is (a, b)
        the interval is between "a (s)" and "b (s)".
        If a is None the beginning of the data is used
        and if b is None then b is set to the end of the interval.
        If baseline is equal to (None, None) all the time
        interval is used. If None, no correction is applied.
    preload : bool
        Preload epochs or not.
    keep_ecg : bool
        When ECG is synthetically created (after picking), should it be added
        to the epochs? Must be False when synthetic channel is not used.
        Defaults to False.
    reject_by_annotation : bool
        Whether to reject based on annotations. If True (default), epochs
        overlapping with segments whose description begins with ``'bad'`` are
        rejected. If False, no rejection based on annotations is performed.

        .. versionadded:: 0.14.0

    verbose : bool, str, int, or None
        If not None, override default verbose level (see :func:`mne.verbose`
        and :ref:`Logging documentation <tut_logging>` for more).

    Returns
    -------
    ecg_epochs : instance of Epochs
        Data epoched around ECG r-peaks.
    """
    has_ecg = 'ecg' in raw or ch_name is not None

    events, _, pulse, ecg = custom_ecg_events(
        raw, custom_args,
        ch_name=ch_name, event_id=event_id, l_freq=l_freq, h_freq=h_freq,
        return_ecg=True, verbose=verbose)

    # Load raw data so that add_channels works
    raw.load_data()

    if not has_ecg:
        ecg_raw = RawArray(
            ecg[None],
            mne.create_info(ch_names=['ECG-SYN'],
                        sfreq=raw.info['sfreq'], ch_types=['ecg']))
        ignore = ['ch_names', 'chs', 'nchan', 'bads']
        for k, v in raw.info.items():
            if k not in ignore:
                ecg_raw.info[k] = v
        raw.add_channels([ecg_raw])

    if keep_ecg:
        if has_ecg:
            raise ValueError('keep_ecg can be True only if the ECG channel is '
                             'created synthetically.')
        else:
            picks = np.append(picks, raw.ch_names.index('ECG-SYN'))
    # create epochs around ECG events and baseline (important)
    ecg_epochs = Epochs(raw, events=events, event_id=event_id,
                        tmin=tmin, tmax=tmax, proj=False, flat=flat,
                        picks=picks, reject=reject, baseline=baseline,
                        reject_by_annotation=reject_by_annotation,
                        verbose=verbose, preload=preload)

    if not has_ecg:
        raw.drop_channels(['ECG-SYN'])

    return ecg_epochs, pulse


@verbose
def custom_bads_ecg(self, inst, custom_args, ch_name=None, threshold=None, start=None,
                  stop=None, l_freq=8, h_freq=16, method='ctps',
                  reject_by_annotation=True, verbose=None):
    """
    Copy-paste from mne.preprocessing.ica, with a few modifications:
        Use qrs_custom() instead of qrs_detector()
        Return pulse
    
    
    Detect ECG related components using correlation.

    .. note:: If no ECG channel is available, routine attempts to create
              an artificial ECG based on cross-channel averaging.

    Parameters
    ----------
    inst : instance of Raw, Epochs or Evoked
        Object to compute sources from.
    ch_name : str
        The name of the channel to use for ECG peak detection.
        The argument is mandatory if the dataset contains no ECG
        channels.
    threshold : float
        The value above which a feature is classified as outlier. If
        method is 'ctps', defaults to 0.25, else defaults to 3.0.
    start : int | float | None
        First sample to include. If float, data will be interpreted as
        time in seconds. If None, data will be used from the first sample.
    stop : int | float | None
        Last sample to not include. If float, data will be interpreted as
        time in seconds. If None, data will be used to the last sample.
    l_freq : float
        Low pass frequency.
    h_freq : float
        High pass frequency.
    method : {'ctps', 'correlation'}
        The method used for detection. If 'ctps', cross-trial phase
        statistics [1] are used to detect ECG related components.
        Thresholding is then based on the significance value of a Kuiper
        statistic.
        If 'correlation', detection is based on Pearson correlation
        between the filtered data and the filtered ECG channel.
        Thresholding is based on iterative z-scoring. The above
        threshold components will be masked and the z-score will
        be recomputed until no supra-threshold component remains.
        Defaults to 'ctps'.
    reject_by_annotation : bool
        If True, data annotated as bad will be omitted. Defaults to True.

        .. versionadded:: 0.14.0

    verbose : bool, str, int, or None
        If not None, override default verbose level (see
        :func:`mne.verbose` and :ref:`Logging documentation <tut_logging>`
        for more). Defaults to self.verbose.

    Returns
    -------
    ecg_idx : list of int
        The indices of ECG related components.
    scores : np.ndarray of float, shape (``n_components_``)
        The correlation scores.

    See Also
    --------
    find_bads_eog

    References
    ----------
    [1] Dammers, J., Schiek, M., Boers, F., Silex, C., Zvyagintsev,
        M., Pietrzyk, U., Mathiak, K., 2008. Integration of amplitude
        and phase statistics for complete artifact removal in independent
        components of neuromagnetic recordings. Biomedical
        Engineering, IEEE Transactions on 55 (10), 2353-2362.
    """
    if verbose is None:
        verbose = self.verbose

    idx_ecg = _get_ecg_channel_index(ch_name, inst)

    if idx_ecg is None:
        if verbose is not None:
            verbose = self.verbose
        ecg, times = _make_ecg(inst, start, stop,
                               reject_by_annotation=reject_by_annotation,
                               verbose=verbose)
    else:
        ecg = inst.ch_names[idx_ecg]

    if method == 'ctps':
        if threshold is None:
            threshold = 0.25
        if isinstance(inst, BaseRaw):
            ecg_epochs, pulse = custom_ecg_epochs(inst, custom_args,
                                                  ch_name=ch_name, keep_ecg=False,
                                                  reject_by_annotation=reject_by_annotation)
            sources = self.get_sources(ecg_epochs).get_data()

            if sources.shape[0] == 0:
                warn('No ECG activity detected. Consider changing '
                     'the input parameters.')
        elif isinstance(inst, BaseEpochs):
            sources = self.get_sources(inst).get_data()
        else:
            raise ValueError('With `ctps` only Raw and Epochs input is '
                             'supported')
        _, p_vals, _ = ctps(sources)
        scores = p_vals.max(-1)
        ecg_idx = np.where(scores >= threshold)[0]
    elif method == 'correlation':
        if threshold is None:
            threshold = 3.0
        scores = self.score_sources(
            inst, target=ecg, score_func='pearsonr', start=start,
            stop=stop, l_freq=l_freq, h_freq=h_freq,
            reject_by_annotation=reject_by_annotation, verbose=verbose)
        ecg_idx = find_outliers(scores, threshold=threshold)
    else:
        raise ValueError('Method "%s" not supported.' % method)
    # sort indices by scores
    ecg_idx = ecg_idx[np.abs(scores[ecg_idx]).argsort()[::-1]]

    self.labels_['ecg'] = list(ecg_idx)
    if ch_name is None:
        ch_name = 'ECG-MAG'
    self.labels_['ecg/%s' % ch_name] = list(ecg_idx)
    return self.labels_['ecg'], scores, pulse
