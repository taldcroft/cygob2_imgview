import Ska.Numpy
import numpy as np
import pywcs
import pyfits

class EventFile(object):
    def __init__(self, filename, hdu=1):
        self.filename= filename
        hdus = pyfits.open(filename)
        self.events = hdus[hdu].data
        self.header = hdus[hdu].header
        self.hdus = hdus

    @property
    def ra_col(self):
        if not hasattr(self, '_ra_col'):
            for key, val in self.header.items():
                if val == 'RA---TAN':
                    self._ra_col = key[5:]
                    break
            else:
                raise ValueError('No RA---TAN ctype found')
        return self._ra_col

    @property
    def dec_col(self):
        if not hasattr(self, '_dec_col'):
            for key, val in self.header.items():
                if val == 'DEC--TAN':
                    self._dec_col = key[5:]
                    break
            else:
                raise ValueError('No DEC--TAN ctype found')
        return self._dec_col

    @property
    def wcs(self):
        """
        Return a pywcs.WCS object corresponding to the given event file HDU.

        :param events_hdu: pyfits HDU
        """
        if not hasattr(self, '_wcs'):
            header = self.header
            ra_col = self.ra_col
            dec_col = self.dec_col
            
            wcs = pywcs.WCS(naxis=2)

            wcs.wcs.equinox = 2000.0
            wcs.wcs.crpix = [header['TCRPX'+ra_col], header['TCRPX'+dec_col]]
            wcs.wcs.cdelt = [header['TCDLT'+ra_col], header['TCDLT'+dec_col]]
            wcs.wcs.cunit = [header['TCUNI'+ra_col], header['TCUNI'+dec_col]]
            wcs.wcs.crval = [header['TCRVL'+ra_col], header['TCRVL'+dec_col]]
            wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

            self._wcs = wcs

        return self._wcs

    def pix2sky(self, i, j):
        return self.wcs.wcs_pix2sky([[i, j]], 1)[0]

    def sky2pix(self, x, y):
        return self.wcs.wcs_sky2pix([[x, y]], 1)[0]

    def event_filter(self, filters):
        if not filters:
            return self.events

        ok = np.ones(len(self.events), dtype=np.bool)
        for colname, limits in filters.items():
            colvals = self.events.field(colname)
            try:
                lo, hi = limits
                if lo is None and hi is not None:
                    ok &= (colvals < hi)
                elif lo is not None and hi is None:
                    ok &= (colvals >= lo)
                elif lo is not None and hi is not None:
                    ok &= (colvals >= lo) & (colvals < hi)
            except TypeError:
                ok &= (colvals == limits)
        return self.events[ok]


    def binned_image(self, x0=None, x1=None, dx=1.0, y0=None, y1=None, dy=1.0, filters=None):
        events = self.event_filter(filters)
        ra_col = self.ra_col
        dec_col = self.dec_col
        dx = float(dx)
        dy = float(dy)

        if x0 is None:
            x0 = np.min(events.field('x'))
        if x1 is None:
            x1 = np.max(events.field('x'))
        if y0 is None:
            y0 = np.min(events.field('y'))
        if y1 is None:
            y1 = np.max(events.field('y'))

        img, x_bins, y_bins = np.histogram2d(events.field('y'), events.field('x'),
                                             bins=[np.arange(y0, y1, dy), np.arange(x0, x1, dx)])

        # Find the position in image coords of the sky pix reference position
        # The -0.5 assumes that image coords refer to the center of the image bin.

        x_crpix = (self.wcs.wcs.crpix[0] - (x0 - dx/2.0)) / dx
        y_crpix = (self.wcs.wcs.crpix[1] - (y0 - dy/2.0)) / dy

        # YIKES, stupid!
        # x_crpix = Ska.Numpy.interpolate(np.arange(len(x_bins)), x_bins, [self.wcs.wcs.crval[0]])[0] - 0.5
        # y_crpix = Ska.Numpy.interpolate(np.arange(len(y_bins)), y_bins, [self.wcs.wcs.crval[1]])[0] - 0.5

        wcs = pywcs.WCS(naxis=2)

        wcs.wcs.equinox = 2000.0
        wcs.wcs.crpix = [x_crpix, y_crpix]
        wcs.wcs.cdelt = [self.wcs.wcs.cdelt[0] * dx, self.wcs.wcs.cdelt[1] * dy, ]
        wcs.wcs.cunit = [self.wcs.wcs.cunit[0], self.wcs.wcs.cunit[1]]
        wcs.wcs.crval = [self.wcs.wcs.crval[0], self.wcs.wcs.crval[1]]
        wcs.wcs.ctype = [self.wcs.wcs.ctype[0], self.wcs.wcs.ctype[1]]

        hdu = pyfits.PrimaryHDU(np.array(img, dtype=np.int32), header=wcs.to_header())
        return hdu

def test(bin=4):
    evt = EventFile('acis_evt2.fits')
    print evt.wcs.wcs_pix2sky([[4096.5, 4096.5]], 1)
    return bin_events(hdus[1], dx=bin, dy=bin)

    
