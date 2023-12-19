import json
import os

import numpy as np
import pandas as pd

import predict_helpers


def predict(FEATURES, MODEL, OUTFILE,  FREEAGO=-7, FEATURES_PASS=None, FREEAGO_PASS=None, KD_CUTOFF=None):

    # check that if passenger strand features are given, a passenger strand free AGO concentration is also given
    if FEATURES_PASS is not None:
        if FREEAGO_PASS is None:
            raise ValueError("Must give free AGO for passenger strand if passenger strand sequences are supplied.")

    # read in features
    features = predict_helpers.process_features(FEATURES, kd_cutoff=KD_CUTOFF)
    features['freeAGO'] = FREEAGO
    if FEATURES_PASS is not None:
        features_pass = predict_helpers.process_features(FEATURES_PASS, kd_cutoff=KD_CUTOFF)
        features_pass['freeAGO'] = FREEAGO_PASS
        features = pd.concat([features, features_pass], sort=False)

    # read in parameters
    with open(MODEL, 'r') as infile:
        TRAIN_PARAMS = json.load(infile)

    FITTED_PARAMS = {'log_decay': TRAIN_PARAMS['log_decay']}
    for feat, val in zip(TRAIN_PARAMS['FEATURE_LIST'][1:], TRAIN_PARAMS['feature_coefs']):
        FITTED_PARAMS[feat + '_coef'] = val

    for param in ['nosite_conc', 'utr3_coef', 'orf_coef']:
        if param in TRAIN_PARAMS:
            FITTED_PARAMS[param] = TRAIN_PARAMS[param]

    FEATURE_LIST = ','.join(TRAIN_PARAMS['FEATURE_LIST'][2:])

    features, predictions = predict_helpers.predict(features, FEATURE_LIST, FITTED_PARAMS)

    features.to_csv(OUTFILE, sep='\t')



