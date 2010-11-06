#!/usr/bin/env python
import sys
import re
import os

import pyfits
import numpy as np
import asciitable
import pyyaks.context
import eventfile
import files_def
import Ska.DBI
import cPickle as pickle

src = pyyaks.context.ContextDict('src')
files = pyyaks.context.ContextDict('files', basedir='/data/cygob2/process/data/baseline/v3')
files.update(files_def.files)

def id_key(cat, id_):
    return '{0}-{1:06d}'.format(cat, int(id_))

def get_detections(filename):
    # lower case, change all whitespace to space, and strip leading/trailing spaces
    lines = [re.sub(r'\s', ' ', x.lower().strip()) for x in open(filename).readlines()]
    dets = asciitable.read(lines, Reader=asciitable.CommentedHeader, numpy=False)
    colnames = dets.dtype.names
    dets_dict = {}
    for det in dets:
        det_dict = dict((name, det[idx]) for idx, name in enumerate(colnames))
        id_ = det_dict['det_id']
        det_dict['cat'] = cat
        det_dict['id'] = id_
        dets_dict[id_key(cat, id_)] = det_dict
    return dets_dict

def get_groups(filename):
    # # Master_number  Group_members  All_matched_flag  Detection_numbers...
    # 1	9	1	1	539	802	6441	6677	6839	6915	17592	17847	
    lines = [re.sub(r'\s', ' ', x.lower().strip()) for x in open(filename).readlines()]
    groups = dict()
    for line in lines[1:]:
        vals = line.split()
        key = int(vals[0])
        groups[key] = {'det_ids': [id_key(cat, x) for x in vals[3:]]}
    return groups
                          
print 'reading'
cat = 'cwdet'
dets = get_detections('all_detections.dat')
groups = get_groups('new_matches.dat')

dbfile = 'datastore.db3'
if os.path.exists(dbfile):
    os.unlink(dbfile)
db = Ska.DBI.DBI(dbi='sqlite', server=dbfile, autocommit=False)
try:
    db.execute('create table dets (key text primary key, val text)')
    db.execute('create table groups (key integer primary key, val text)')
except:
    pass

print 'writing dets'
for key in sorted(dets):
    db.insert({'key': key, 'val': pickle.dumps(dets[key])}, 'dets')
db.commit()
        
print 'writing groups'
for key in sorted(groups):
    db.insert({'key': key, 'val': pickle.dumps(groups[key])}, 'groups')
db.commit()


