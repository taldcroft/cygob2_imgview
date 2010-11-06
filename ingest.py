#!/usr/bin/env python
import os

import vots
import Ska.DBI
import cPickle as pickle

def id_key(cat, id_):
    return '{0}-{1}'.format(cat, id_)

dbfile = 'datastore.db3'
if os.path.exists(dbfile):
    os.unlink(dbfile)
db = Ska.DBI.DBI(dbi='sqlite', server=dbfile, autocommit=False)
for table_name, id_type in (('dets', 'text'),
                            ('groups', 'int')):
    db.execute("""CREATE TABLE {0} (
                  id {1} PRIMARY KEY,
                  ra float,
                  dec float,
                  flag int,
                  info text)""".format(table_name, id_type))
print 'reading groups'
hdr, groups = vots.parse_table(open('cat_groups.vots'))
        
print 'reading assoc_groups_obsids'
hdr, groups_obsids = vots.parse_table(open('assoc_groups_obsids.vots'))
groups_obsids_dict = {}
for group_id, obsid in groups_obsids:
    group_id = int(group_id)
    obsid = int(obsid)
    groups_obsids_dict.setdefault(group_id, []).append(obsid)

print 'reading assoc_groups'
cats = ('cwdet', 'pwdet', 'vdet', 'mopwdet', 'hdet')
hdr, groups_detids = vots.parse_table(open('assoc_groups.vots'))
groups_detids_dict = {}
for group_id, cat_id, detid, flag, conf in groups_detids:
    group_id = int(group_id)
    detid = int(detid)
    groups_detids_dict.setdefault(group_id, []).append(id_key(cats[cat_id], detid))

for cat in ('cwdet', 'pwdet'):
    print 'reading {0} dets'.format(cat)
    hdr, dets = vots.parse_table(open('cat_{0}.vots'.format(cat)))

    print 'writing {0} dets'.format(cat)
    for det in dets:
        info = dict((x, det[x].tolist()) for x in det.dtype.names)
        info['id'] = id_key(cat, det['id'])
        db.insert({'id': info['id'],
                   'ra': det['ra'],
                   'dec': det['dec'],
                   'flag': 0,
                   'info': pickle.dumps(info)}, 'dets')
    db.commit()

print 'writing groups'
for i, group in enumerate(groups):
    try:
        info = dict((x, group[x].tolist()) for x in group.dtype.names)
        info['obsids'] = groups_obsids_dict[group['id']]
        info['det_ids'] = groups_detids_dict[group['id']]
        info['views'] = []
        info['comment'] = ""
        db.insert({'id': group['id'],
                   'ra': group['ra'],
                   'dec': group['dec'],
                   'flag': 0,
                   'info': pickle.dumps(info)}, 'groups')
    except Exception, err:
        print 'Failure for group_id={0}: {1}'.format(group['id'], err)
        pass
db.commit()

for table_name, id_type in (('dets', 'text'),
                            ('groups', 'int')):
    for col_name in ('ra', 'dec'):
        db.execute("CREATE INDEX idx_{0}_{1} ON {0} ( {1} )".format(
            table_name, col_name))

db.commit()
