#!/usr/bin/env python

import synapseclient
import os
from os.path import join

target_dir = '/Users/daiwei89/storage/cs_proj/dap/data/KIRC'

# GBM
#synapse_ids = ['syn1710366', 'syn1710372', 'syn1710374', 'syn1710368', 'syn1710370', 'syn1715822']

# KIRC
synapse_ids = ['syn1714090', 'syn1714093', 'syn1710287', 'syn1710306', 'syn1710293',
    'syn1710289', 'syn1710291', 'syn1710303', 'syn1715824']

def MoveFile(entity):
    assert len(entity['files']) == 1
    src = join(entity['cacheDir'], entity['files'][0])
    cmd = 'mv %s %s' % (src, target_dir)
    print cmd
    os.system(cmd)

syn = synapseclient.Synapse()
syn.login(email='daiwei89@gmail.com', password='david123')

os.system('mkdir -p %s' % target_dir)

for synapse_id in synapse_ids:
  print 'Processing', synapse_id
  MoveFile(syn.get(synapse_id))
