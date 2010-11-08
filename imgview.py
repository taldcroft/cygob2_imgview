#!/usr/bin/env python
import sys
import os
import time
import re
import pygtk
import gtk
import cPickle as pickle

import Ska.DBI
import pyfits
import aplpy
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import asciitable
import pyyaks.context
import eventfile

import files_def

src = pyyaks.context.ContextDict('src')
files = pyyaks.context.ContextDict('files', basedir='/data/cygob2/process/data/baseline/v3')
files.update(files_def.files)
dbfile = '/data/cygob2/tables/datastore.db3'

geom = dict(main_image_size=8,
            main_image_arcsec=30,
            n_obs_frames=2)
bands = dict(broad=(500, 7000),
             soft=(500, 1200),
             medium=(1200, 2000),
             hard=(2000,7000))

class ImageDisplay(object):
    def __init__(self, win, size, show_labels=True):
        fig = Figure(figsize=(size, size))
        canvas = FigureCanvas(fig)  # a gtk.DrawingArea
        canvas.set_size_request(500, 500)
        toolbar = NavigationToolbar(canvas, win)
        self.vbox = gtk.VBox()
        self.vbox.pack_start(canvas)
        self.vbox.pack_start(toolbar, False, False)
        self.fig = fig
        self.show_labels = show_labels
        self.eventfiles = {}

    def refresh(self):
        self.gc.refresh()

    def update_image(self, image):
        self.fig.clear()
        gc = aplpy.FITSFigure(image, figure=self.fig, auto_refresh=False)
        gc.show_colorscale(cmap='hot', vmin=0, vmax=np.max(image.data))
        gc.add_colorbar()
        gc.axis_labels.hide()
        gc.tick_labels.set_xformat('ddd.dddd')
        gc.tick_labels.set_yformat('ddd.dddd')
        gc.tick_labels.set_font(size='small')
        self.gc = gc

    def update_regions(self, regions):
        for i, region in enumerate(regions):
            linewidth = (4.0 if region['active'] else 1.0)
            self.gc.show_circles([region['ra']], [region['dec']], [region['radius']],
                                 edgecolor=region['edgecolor'], linewidth=linewidth,
                                 layer='circle_set_%d' % (i+1))


class TableColumn(dict):
    def __init__(self, table, name, col):
        self.table = table
        self.name = name
        self.col = col
        self.table.attach(gtk.Label(name.title()), self.col, self.col + 1, 0, 1)
        dict.__init__(self)

    def __setitem__(self, row, widget):
        if row in self:
            self.table.remove(self[row])
        dict.__setitem__(self, row, widget)
        self.table.attach(widget, self.col, self.col + 1, row + 1, row + 2, xoptions=gtk.SHRINK)
        widget.show()

    def clear(self):
        for widget in self.values():
            self.table.remove(widget)
        dict.clear(self)


class DetsTable(object):
    def __init__(self):
        self.colnames = 'id obsid band ra dec psf_size net_counts src_significance'.split()
        self.fmt = dict(ra='{0:.4f}', dec='{0:.4f}', psf_size='{0:.2f}')
        self.n_rows = 1
        self.n_cols = 1 + len(self.colnames)

        self.table = gtk.Table(rows=self.n_rows, columns=self.n_cols)

        self.vbox = gtk.VBox(False, 0)
        self.vbox.pack_start(self.table, False, False, 0)

        self.scrolled_window = gtk.ScrolledWindow()
        self.scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.scrolled_window.add_with_viewport(self.vbox)

        self.region_select = TableColumn(self.table, 'region', 0)
        self.det_labels = {}
        for j, colname in enumerate(self.colnames):
            self.det_labels[colname] = TableColumn(self.table, colname, 1 + j)
            
    def update(self, group_dets):
        self.region_select.clear()
        for det_labels in self.det_labels.values():
            det_labels.clear()

        self.n_rows = len(group_dets) + 1
        self.table.resize(self.n_rows, self.n_cols)
        
        for i, det in enumerate(group_dets):
            for colname in self.colnames:
                strval = self.fmt.get(colname, '{0}').format(det[colname])
                self.det_labels[colname][i] = gtk.Label(strval)
            self.region_select[i] = gtk.CheckButton()


class ValsRadio(object):
    def __init__(self, name):
        self.widgets = []
        self.hbox = gtk.HBox(False, 0)
        self.name = name
        self.vals = []

    def update(self, vals, val_select=None):
        self.vals = list(vals)
        for widget in self.widgets:
            self.hbox.remove(widget)
        self.widgets = []

        for j, val in enumerate(self.vals):
            group = (self.widgets[0] if j else None)
            widget = gtk.RadioButton(group=group, label=str(val).title())
            self.widgets.append(widget)
            self.hbox.pack_start(widget, False, False, 0)
            if val_select is not None:
                widget.set_active(val == val_select)
            widget.show()

        for widget in self.widgets:
            widget.connect('toggled', self.callback)

    def get_val(self):
        for val, widget in zip(self.vals, self.widgets):
            if widget.get_active():
                return val
        else:
            raise ValueError('No radio button selected for {0}'.format(self.name))

    def set_callback(self, callback):
        self.callback = callback

class InfoPanel(object):
    def __init__(self):
        self.vbox = gtk.VBox(False, 0)
        # Top row of nav buttons
        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_GO_FORWARD, gtk.ICON_SIZE_BUTTON)
        self.next = gtk.Button()
        self.next.set_image(image)

        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_GO_BACK, gtk.ICON_SIZE_BUTTON)
        self.prev = gtk.Button()
        self.prev.set_image(image)

        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_APPLY, gtk.ICON_SIZE_BUTTON)
        self.ok_button = gtk.Button(label='ok')
        self.ok_button.set_image(image)

        image = gtk.Image()
        image.set_from_stock(gtk.STOCK_CANCEL, gtk.ICON_SIZE_BUTTON)
        self.bad_button = gtk.Button(label='bad')
        self.bad_button.set_image(image)

        self.randomize = gtk.ToggleButton('Randomize')
        self.randomize.set_active(False)
        self.group_id_entry = gtk.Entry(max=6)
        self.group_id_entry.set_width_chars(6)


        hbox1 = gtk.HBox(False, 0)
        hbox1.pack_start(self.prev, False, False, 0)
        hbox1.pack_start(self.next, False, False, 0)
        hbox1.pack_start(self.ok_button, False, False, 0)
        hbox1.pack_start(self.bad_button, False, False, 0)
        hbox1.pack_start(self.group_id_entry, False, False, 0)
        hbox1.pack_start(self.randomize, False, False, 0)

        hbox2 = gtk.HBox(False, 0)
        self.obsids_radio = ValsRadio('obsids')
        hbox2.pack_start(self.obsids_radio.hbox, False, False, 0)
        
        hbox3 = gtk.HBox(False, 0)
        self.bands_radio = ValsRadio('bands')
        hbox3.pack_start(self.bands_radio.hbox, False, False, 0)
        
        hbox4 = gtk.HBox()
        hbox4.set_border_width(4)
        sw = gtk.ScrolledWindow()
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.comment_textview = gtk.TextView()
        self.comment_textview.set_wrap_mode(gtk.WRAP_WORD)
        self.comment_textbuffer = self.comment_textview.get_buffer()
        sw.add(self.comment_textview)
        hbox4.pack_start(sw, True, True, 3)

        sw = gtk.ScrolledWindow()
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        self.views_textview = gtk.TextView()
        self.views_textview.set_wrap_mode(gtk.WRAP_WORD)
        self.views_textbuffer = self.views_textview.get_buffer()
        sw.add(self.views_textview)
        hbox4.pack_start(sw, True, True, 3)

        self.dets_table = DetsTable()
        self.other_dets_table = DetsTable()
        self.vbox.pack_start(hbox1, False, False, 0)
        self.vbox.pack_start(hbox2, False, False, 0)
        self.vbox.pack_start(hbox3, False, False, 0)
        self.vbox.pack_start(hbox4, False, False, 0)
        self.vbox.pack_start(self.dets_table.scrolled_window, True, True, 0)
        self.vbox.pack_start(self.other_dets_table.scrolled_window, True, True, 0)

    def get_group_id(self):
        return int(self.group_id_entry.get_text())

    def set_group_id(self, group_id):
        self.group_id_entry.set_text(str(group_id))

    group_id = property(get_group_id, set_group_id)


class Controller(object):
    def __init__(self, info_panel, image_display):
        self.db = Ska.DBI.DBI(server=dbfile, dbi='sqlite', autocommit=False, numpy=False)
        self.dets, self.dets_table = self.fetch_detections()
        self.det_ras = np.array([x['ra'] for x in self.dets_table])
        self.det_decs = np.array([x['dec'] for x in self.dets_table])
        self.groups = self.fetch_groups()
        self.group_ids = [x['id'] for x in self.db.fetchall('select id from groups order by id asc')]
        self.group = None
        self.index = 0
        self.image_display = image_display
        self.info_panel = info_panel
        self.evt2_cache = {}
        self.image_cache = {}

    def new_group(self, widget=None, offset=None, status=None):
        # Store comment (if modified) to database for current group
        self.store_group(status)
        
        if offset is None:
            try:
                self.index = self.group_ids.index(self.info_panel.group_id)
            except ValueError:
                pass
        else:
            self.index = (self.index + offset) % len(self.group_ids)

        group_id = self.group_ids[self.index]
        self.info_panel.group_id = group_id
        self.group = self.fetch_group(group_id)
        self.groups[group_id] = self.group  # refresh the in-memory version
        det_ids = sorted(self.group['info']['det_ids'],
                         key=lambda x: -self.dets[x]['src_significance'])
        self.group_dets = [self.dets[x] for x in det_ids]

        # Find detections near image window
        ra0 = self.group['ra']
        dec0 = self.group['dec']
        dec_halfwidth = geom['main_image_arcsec'] / 3600. / 2.
        ra_halfwidth = dec_halfwidth / np.cos(np.radians(dec0))
        ok = ((abs(self.det_ras - ra0) < ra_halfwidth)
              & (abs(self.det_decs - dec0) < dec_halfwidth))
        self.other_dets = [self.dets_table[i] for i in np.flatnonzero(ok)
                           if self.dets_table[i]['id'] not in det_ids]

        views_text = '\n'.join('{0:15s} {1:15s} {2}'.format(x['date'], x['user'],
                                                              x['status'] or '-')
                               for x in self.get_user_views())
        self.info_panel.views_textbuffer.set_text(views_text)
        self.info_panel.comment_textbuffer.set_text(self.group['info']['comment'])

        obsids = sorted(self.group['info']['obsids'])
        self.info_panel.obsids_radio.update(obsids, self.dets[det_ids[0]]['obsid'])
        self.info_panel.bands_radio.update(['broad', 'soft', 'medium', 'hard'])

        self.info_panel.dets_table.update(self.group_dets)
        for widget in info_panel.dets_table.region_select.values():
            widget.connect('toggled', self.update_regions)

        self.info_panel.other_dets_table.update(self.other_dets)
        for widget in info_panel.other_dets_table.region_select.values():
            widget.connect('toggled', self.update_regions)

        self.update_image(None)

    def update_image(self, widget):
        # RadioButton group toggle fires two callbacks, just respond to active one
        if isinstance(widget, gtk.RadioButton) and not widget.get_active():
            return 
        ra, dec = self.group['ra'], self.group['dec']
        obsid = self.info_panel.obsids_radio.get_val()
        band = self.info_panel.bands_radio.get_val()
        if src['obsid'].val != obsid:
            src['obsid'] = obsid
            filename = files['evt2.fits'].abs
            if filename not in self.evt2_cache:
                print 'Reading event file', filename
                self.evt2_cache[filename] = eventfile.EventFile(filename)
            self.evt2 = self.evt2_cache[filename]
        x, y = self.evt2.sky2pix(ra, dec)
        try:
            key = (obsid, band, x, y)
            image = self.image_cache[key]
        except KeyError:
            sz = geom['main_image_arcsec'] # divide by 2 for half-width and mult by 2 for pixels
            image = self.evt2.binned_image(x0=x-sz, x1=x+sz, y0=y-sz, y1=y+sz,
                                           filters={'energy': bands[band]})
            self.image_cache[key] = image
        self.image_display.update_image(image)
        self.update_regions()

    def update_regions(self, widget=None):
        regions = []
        for dets, dets_table, edgecolor in (
            (self.group_dets, self.info_panel.dets_table, 'g'),
            (self.other_dets, self.info_panel.other_dets_table, 'm')):
            for i, det in enumerate(dets):
                regions.append(dict(ra=det['ra'],
                                    dec=det['dec'],
                                    radius=det['psf_size'] * 0.5 / 3600.0,
                                    active=dets_table.region_select[i].get_active(),
                                    edgecolor=edgecolor))
        self.image_display.update_regions(regions)
        self.image_display.refresh()

    def update_randomize(self, widget):
        if widget.get_active():
            np.random.shuffle(self.group_ids)
        else:
            self.group_ids = sorted(groups.keys())
        self.index = self.group_ids.index(self.info_panel.group_id)

    def fetch_detections(self):
        # lower case, change all whitespace to space, and strip leading/trailing spaces
        dets = self.db.fetchall('select * from dets')
        dets_dict = {}
        for det in dets:
            det['info'] = pickle.loads(str(det['info']))
            det.update(det['info'])
            dets_dict[det['id']] = det
        return dets_dict, dets

    def fetch_groups(self):
        groups = self.db.fetchall('select * from groups')
        groups_dict = {}
        for group in groups:
            group['info'] = pickle.loads(str(group['info']))
            groups_dict[group['id']] = group
        return groups_dict

    def fetch_group(self, group_id):
        group = self.db.fetchone('select * from groups where id={0}'.format(group_id))
        if group is None:
            raise ValueError('No group matching id={0}'.format(group_id))
        group['info'] = pickle.loads(str(group['info']))
        return group

    def get_user_views(self):
        user_views = {}
        for view in self.group['info']['views']:
            user = view['user']
            if user not in user_views or view['status'] is not None:
                user_views[user] = view
        return sorted(user_views.values(), key=lambda x: -x['time'])

    def store_group(self, status=None):
        if self.group is None:
            return

        buff = self.info_panel.comment_textbuffer
        start, end = buff.get_bounds()
        self.group['info']['comment'] = buff.get_text(start, end)

        # Add the current view and status info
        newview = {'user': os.environ['USER'],
                   'time': time.time(),
                   'date': time.strftime('%Y-%m-%d %H:%M'),
                   'status': status}
        self.group['info']['views'].append(newview)

        group = self.group.copy()
        group['info'] = pickle.dumps(group['info'])
        self.db.insert(group, 'groups', replace=True, commit=True)


# Create the main window
win = gtk.Window()
win.connect("destroy", lambda x: gtk.main_quit())
win.set_default_size(1200, 700)
win.set_title("CygOB2 Image Browser")

image_display = ImageDisplay(win, size=geom['main_image_size'])
info_panel = InfoPanel()

main_box = gtk.HBox(homogeneous=False, spacing=0)
main_box.pack_start(info_panel.vbox)
main_box.pack_start(image_display.vbox, False, False, 0)
win.add(main_box)

controller = Controller(info_panel, image_display)

# Set up signals
info_panel.next.connect('clicked', controller.new_group, 1)
info_panel.prev.connect('clicked', controller.new_group, -1)
info_panel.ok_button.connect('clicked', controller.new_group, 1, 'OK')
info_panel.bad_button.connect('clicked', controller.new_group, 1, 'BAD')
info_panel.group_id_entry.connect('activate', controller.new_group)
info_panel.randomize.connect('toggled', controller.update_randomize)
info_panel.obsids_radio.set_callback(controller.update_image)
info_panel.bands_radio.set_callback(controller.update_image)

controller.new_group()
win.show_all()
gtk.main()
