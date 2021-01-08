import cfgrib
import rasterio
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import sys

class HrrrSeries:
    def __init__(self, fnames, fo, i0, j0, i1, j1):
        self.fnames = fnames
        self.fo = fo
        self.i0 = i0
        self.j0 = j0
        self.i1 = i1
        self.j1 = j1

        # load the first data
        hr = Hrrr(fnames[0], fo, i0, j0, i1, j1)
        self.hr0 = hr
        self.ymdh0 = hr.ymdh
        self.ni = hr.ni
        self.nj = hr.nj
        self.nhrs = len(fnames) # assume hourly
        self.xlongdot = hr.xlongdot
        self.xlatdot = hr.xlatdot

        self.crs = hr.crs

        self.gtrans0 = hr.gtrans
        self.dx = hr.dx
        self.ni0 = hr.ni0
        self.nj0 = hr.nj0
        self.nk0 = hr.nk0



    def write_header(self):
        self.hr0.write_header(self.nhrs)

    def write_data(self):
        self.hr0.write_data()

        for fname in self.fnames[1:]:
            hr = Hrrr(fname, self.fo, self.i0, self.j0, self.i1, self.j1)





class Hrrr:

    def __init__(self, fname, fo, i0, j0, i1, j1):
        """
        read grib2 file for hrrr and writes 3d.dat file 

        :param fname: input, subsetted, grib file
        :param fo: output file handle
        :param i0: 1-based easting index of ll coner grid cell of extracted region
        :param j0: 1-based norghing index of ll coner grid cell of extracted region
        :param i1: 1-based easting index of ur coner grid cell of extracted region
        :param j1: 1-based norghing index of ur coner grid cell of extracted region
        """
        
        # Ouptut text file
        try:
            if not fo.writable:
                raise RuntimeError("output file not writable")
        except AttributeError:
            fo = open(fo, 'w')
        self.fo = fo

        # LL corner index, 1-based
        self.i0=i0
        self.j0=j0


        # Open un-extracted grib file for spatial info
        r = rasterio.open(fname)

        # Read spatial info for the original grib file
        self.shape0 = r.shape
        self.nj0, self.ni0 = r.shape
        self.crs = r.crs.to_dict()
        self.gtrans = r.get_transform()
        assert self.gtrans[2] == 0 and self.gtrans[4] == 0
        self.x0 = self.gtrans[0]
        self.y0 = self.gtrans[3]
        if self.gtrans[5] < 0: 
            self.y0 += self.gtrans[5] * (self.nj0 - 1)
        assert abs(self.gtrans[1]) == abs(self.gtrans[5])
        self.dx = self.gtrans[1]

        # subset grib data, create temporary file
        #fname2 = 'tmp.grib2'
        fid, fname2 = tempfile.mkstemp()
        self.fname2 = fname2
        ext_grib(fname, fname2, (i0,j0,i1,j1))

        # Open data for reading
        ds = cfgrib.open_datasets(fname2)

        # Read common data

        # time stamp
        self.ymdh = pd.to_datetime(ds[3].time.values).strftime('%Y%m%d%H')

        # sigma levels
        self.sigma = ds[2].isobaricInhPa.values / 1013
        self.nk0 = len(self.sigma)

        # horizontal index
        self.nj, self.ni = ds[3].sp.shape
        jdx,idx = np.indices(ds[3].sp.shape)
        self.jdx = jdx
        self.idx = idx
        self.iindex = idx + i0
        self.jindex = jdx + j0

        # horozonta coords
        self.xlatdot = ds[3].latitude
        self.xlongdot = ds[3].longitude - 360

        # elevation
        self.ielevdot = ds[3].orog

        # mm5 land use
        #self.iland = ds[3].lsm
        self.iland = -9 * np.ones_like(ds[3].lsm)

        # Read surface data

        # surface pressure, Pa => hPa
        self.spres = ds[3].sp / 100 

        # precipitation, kg m-2 => cm
        self.rain = ds[3].tp / 10 

        # snow cover, 1 for >0, 0 fo 0
        self.sc = ds[3].snowc  > 0

        # short and long wave radiation
        self.radsw = ds[3].dswrf
        self.radlw = ds[3].dlwrf

        # T 2m
        self.t2 = ds[1].t2m

        # specific hum 2m kg/kg => g/kg
        self.q2 = ds[1].sh2 / 1000

        # wd ws
        u10 = ds[0].u10
        v10 = ds[0].v10
        self.ws10, self.wd10 = uv2ws(u10, v10)

        # surface T
        self.sst = ds[3].t

        # Read upper data

        # pressure
        self.pres = ds[2].isobaricInhPa.values

        # hiehgt
        self.z = ds[2].gh

        # temp
        self.tempk = ds[2].t

        # ws wd
        u = ds[2].u
        v = ds[2].v
        self.ws, self.wd = uv2ws(u, v)

        # relative humidity (%)
        self.rh = ds[2].r

        # vertical velocity, vaper mix ratio
        w = ds[2].w # w in Pa/sec
        shp = w.shape
        self.w, self.vapmr = getw(w, np.repeat(self.pres,
            np.prod(shp[1:])).reshape(shp), self.tempk, self.rh)


        # cloud mixing ratio
        self.cldmr = ds[2].clwmr

        #print('w', ds[2].w.units)
        #print('r', ds[2].r.units)
        #print('gh', ds[2].gh.units)
        #print('t', ds[2].t.units)
        #print('clwmr', ds[2].clwmr.units)
        #
        #print('pres', [_(self.pres) for _ in (np.min, np.max, np.mean)])
        #print('w', [_(w.values) for _ in (np.min, np.max, np.mean)])
        #print('tempk', [_(self.tempk.values) for _ in (np.min, np.max, np.mean)])
        #print('rh', [_(self.rh.values) for _ in (np.min, np.max, np.mean)])
        #print('w', [_(self.w.values) for _ in (np.min, np.max, np.mean)])
        #print('cldmr', [_(self.cldmr.values) for _ in (np.min, np.max, np.mean)])

        ## all necessary data got read.  clean up
        #del ds
        #os.unlink(fname2)


    def write_header(self, nhrs=1):

        # header 1
        self.fo.write( '%-132s\n' % 
                '3D.DAT          2.1             Header Structure with Comment Lines'
                )
        # header 2
        self.fo.write('   1\n') 
        self.fo.write('%-132s\n' % 
                'Produced by hrrr2m3d v0.8'
                )

        # header 3
        self.fo.write( '  1  1  1  0  0  0' + '\n')
        
        # header 4
        # use rasterio
        assert self.crs['proj'] == 'lcc'
        template = '%-4s%9.4f%10.4f%7.2f%7.2f%10.3f%10.3f%8.3f%4d%3d%3d\n'
        self.fo.write(
                template % (
                    self.crs['proj'].upper(),
                    self.crs['lat_0'],
                    self.crs['lon_0'] - 360,
                    self.crs['lat_1'],
                    self.crs['lat_2'],
                    self.x0 / 1000,
                    self.y0 / 1000,
                    self.dx / 1000,
                    self.ni0,
                    self.nj0,
                    self.nk0,
                    )
                )


        # header 5
        self.fo.write( 
                '  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0'
                + '\n'
                )

        # header 6
        template = '%10s%5d%4d%4d%4d\n'
        self.fo.write(
                template % (
                    self.ymdh,
                    nhrs,
                    self.ni,
                    self.nj,
                    len(self.pres),
                    )
                )

        # header 7
        template = '%4d%4d%4d%4d%4d%4d%10.4f%10.4f%9.4f%9.4f\n'
        self.fo.write( 
                template % (
                    self.i0,
                    self.j0,
                    self.i0 + self.ni - 1,
                    self.j0 + self.nj - 1,
                    1,
                    len(self.pres), # TODO check if calmet care
                    self.xlongdot.min(), # TODO, is dot ok?
                    self.xlongdot.max(),
                    self.xlatdot.min(),
                    self.xlatdot.max(),
                    )
                )

        # header for upper layers
        lines = ['%6.3f' % _ for _ in self.sigma]
        self.fo.write('\n'.join(lines) + '\n')

        # header for surface
        for ii,jj in zip(self.idx.flat, self.jdx.flat):
            line = ''
            line += '%4d%4d' % (self.iindex[jj,ii], self.jindex[jj,ii])
            line += '%9.4f' % self.xlatdot[jj,ii]
            line += '%10.4f' % self.xlongdot[jj,ii]
            line += '%5d' % self.ielevdot[jj,ii]
            line += '%3d' % self.iland[jj,ii]
            line += ' %9.4f%10.4f%5d' % (-999,-999,-999)
            self.fo.write(line + '\n')


    def write_data(self):
        for ii,jj in zip(self.idx.flat, self.jdx.flat):
            self.write_ij(ii, jj)


    def write_ij(self, ii, jj):
        self.write_surface( ii,jj)
        self.write_upper( ii,jj)


    def write_upper(self,  ii, jj):
        template = '%4d%6d%6.1f%4d%5.1f%6.2f%3d%5.2f%6.3f' 

        for k in range(len(self.pres)):

            self.fo.write(template % ( 
                self.pres[k], 
                self.z[k,jj,ii],
                self.tempk[k,jj,ii],
                self.wd[k,jj,ii],
                self.ws[k,jj,ii],
                self.w[k,jj,ii], 
                self.rh[k,jj,ii],
                self.vapmr[k,jj,ii], 
                self.cldmr[k,jj,ii],
                ) + '\n')


    def write_surface(self,  ii, jj):
        template = '%10s%03d%03d%7.1f%5.2f%2d%8.1f%8.1f%8.1f%8.2f%8.1f%8.1f%8.1f'
        line = template % (
            self.ymdh, 
            self.iindex[jj,ii], 
            self.jindex[jj,ii], 
            self.spres[jj,ii], 
            self.rain[jj,ii], 
            self.sc[jj,ii], 
            self.radsw[jj,ii], 
            self.radlw[jj,ii], 
            self.t2[jj,ii], 
            self.q2[jj,ii], 
            self.wd10[jj,ii], 
            self.ws10[jj,ii], 
            self.sst[jj,ii],
                )
        self.fo.write(line + '\n')


def esat(tdegc):
    """
    compute the saturation water vapor pressure (mb)

    using the method of Lowe (1977) (JAM, 16, pp 100-103).

    :params tdegc: air temperature (deg C)

    :returns esat: saturation water vapor pressure (mb)
    """
    a = np.array([6.107799961,
            4.436518521e-1,
            1.428945805e-2, 
            2.650648471e-4, 
            3.031240396e-6, 
            2.034080948e-8, 
            6.136820929e-11])
    esat = a[0] + tdegc * (
            a[1] + tdegc * (
                a[2] + tdegc *( 
                    a[3] + tdegc *( 
                        a[4] + tdegc *( 
                            a[5] + tdegc *a[6]
                            )))))
    esat = np.maximum(0, esat)
    return esat

def gettv(p, tk, rh):
    """
    get virtual temperature

    :param p: pressure in mb
    :param tk: air temperature in K
    :param rh: relative humidity in %

    :returns:
        - tv - virtual temperature in K
        - aq - vapor mixing ratio in g/kg
    """
    

    #assert all(rh <= 100)
    rh = np.minimum(rh, 100)

    tc = tk - 273.15
    es0 = esat(tc)
    ee = es0 * rh / 100

    qq = .622 * ee / (p - .378 * ee)
    qq = np.maximum(qq, 0)

    tv = tk * (1 + .608 * qq)
    aq = qq * 1000
    return tv, aq


def getw(wp, pp, atk, rh):
    """
    convert omega velocity to w velocity

    # TODO doc taken from here, check accuracy!!
    # https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml
    :param wp: omega in Pa/s
    :param pp: pressure in hPa
    :param atk: temperature in K
    :param rh: relative humidity

    :returns:
        - ww - w in m/s
        - aq - vapor mixing ratio in g/kg
    """
    g = 9.80616
    Rd = 287.0
    #alpha = Rd/g

    tv, aq = gettv(pp, atk, rh)

    ro = pp * 100 / Rd / tv
    ww = -wp / ro / g

    return ww, aq


def uv2ws(u,v):
    """
    wind speed/direction from u,v

    :param u: u component (m/sec)
    :param v: v component (m/sec)

    :returns:
        - ws - wind speed in m/s
        - wd - wind direction (0 deg = due north, clockwise)
    """
    ws = np.sqrt(u*u + v*v)
    wd = 270. - 180/np.pi * np.arctan2(v, u)
    return ws,wd


def ext_grib(fname, fname2, rng):
    """
    make smaller grib by subsetting variables and horizontal extent

    :param fname: input grib file name
    :param rng: horitontal extent, 1based index of llidx, lljdx, uridx, urjdx

    """

    # horizontal extent
    i0, j0, i1, j1 = rng

    # variables of interest
    vnames = ('PRES', 'APCP', 'SNOWC', 'DSWRF', 'DLWRF', 'TMP', 'SPFH', 'UGRD',
    'VGRD', 'WIND', 'HGT', 'VVEL', 'RH', 'CLMR', 'LAND')

    # subset by location
    irng = f'{i0}:{i1}'
    jrng = f'{j0}:{j1}'

    # subset by variables
    # all levels of intrest
    levels = ['.* mb' ]
    levels += ['surface']
    levels += [f"{z} m above ground" for z in [2,10]]
    lev = ':(' + '|'.join(levels) + '):'

    # all variables of interest
    nam = ':(' + '|'.join(vnames) + '):'

    cmd = ['wgrib2', '-v0',
            fname, 
            '-if', lev, 
            '-if', nam, 
            '-ijsmall_grib', irng, jrng,
            fname2]
     
    print(cmd)
    subprocess.run(cmd)

def tester():
    # test file
    fname = 'hrrr.t00z.wrfprsf00.grib2'

    # ector_midland
    i0, i1 = 730, 763
    j0, j1 = 277, 295

    hrrr = Hrrr(fname, 'test.txt', i0, j0,i1, j1)


    hrrr.write_header()

    hrrr.write_data()


if __name__ == '__main__':
    import sys
    fname, oname, i0, j0, i1, j1 = sys.argv[1:7]
    i0, j0, i1, j1 = [int(_) for _ in (i0, j0, i1, j1)]

    hrrr = Hrrr(fname, oname, i0, j0, i1, j1)
    hrrr.write_header()
    hrrr.write_data()


