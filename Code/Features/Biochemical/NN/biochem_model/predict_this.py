#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json
from optparse import OptionParser
import os

import numpy as np
import pandas as pd

import predict_helpers


def predict_this(FEATURES,MODEL,FREEAGO,OUTFILE):
    KD_CUTOFF=0.0
    
    # read in features
    features = predict_helpers.process_features(FEATURES, kd_cutoff=KD_CUTOFF)
    features['freeAGO'] = FREEAGO

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

    _, predictions = predict_helpers.predict(features, FEATURE_LIST, FITTED_PARAMS)

    #predictions.to_csv(OUTFILE, sep='\t')
    return predictions



