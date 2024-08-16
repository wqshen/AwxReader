# -*- coding: utf-8 -*-
# @Author: wqshen
# @Email: wqshen91@gmail.com
# @Date: 2023/2/14 16:22
# @Last Modified by: wqshen


import os
import sys
import argparse
import numpy as np
import xarray as xr
from os import PathLike
from datetime import datetime
from struct import calcsize, unpack
from pyproj import CRS, Transformer
from typing import cast, Optional, Union
from dataclasses import dataclass, fields
from .awx_geos_image import AwxGeosImageHead
from .awx_polar_image import AwxPolarImageHead
from .awx_grid_field import AwxGridHead
from .awx_discrete_field import AwxDiscreteHead


@dataclass
class AwxFileName:
    satellite: str
    product_type: str
    channel: str
    proj: str
    obs_date: str
    obs_start_time: str

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            setattr(self, field.name, field.type(field_value))


@dataclass
class AwxDaasFileName:
    z: str
    dtype: str
    c: str
    center: str
    gentime: str
    p: str
    satellite: str
    product_type: str
    channel: str
    proj: str
    obs_date: str
    obs_start_time: str

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            setattr(self, field.name, field.type(field_value))


@dataclass
class AwxHead:
    filename_Sat96: str
    int_order: int = 0
    head1_length: int = 40
    head2_length: int = 0
    filled_length: int = 0
    record_length: int = 0
    head_record: int = 40
    data_record: int = 0
    product_kind: int = 0
    zip_method: int = 0
    form_string: str = "SAT2004"
    data_quality: int = 0

    @staticmethod
    def dtype() -> str:
        return '<12s9h8sh'

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            if isinstance(field_value, bytes):
                field_value = field_value.rstrip(b'\x00').decode('utf8')
            setattr(self, field.name, field.type(field_value))


@dataclass
class AwxPositioning:
    coordination: int = 0
    source: int = 0
    resolution: int = 0
    ul_lat: int = 0
    ul_lon: int = 0
    size_h: int = 0
    size_v: int = 0
    reserve: int = 0

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            if isinstance(field_value, bytes):
                field_value = field_value.rstrip(b'\x00').decode('utf8')
            setattr(self, field.name, field.type(field_value))


@dataclass
class Awx(object):
    pathfile: Optional[Union[str, PathLike]] = None
    name: Optional[AwxFileName | AwxDaasFileName] = None
    head1: AwxHead = None
    head2: Union[AwxGeosImageHead, AwxPolarImageHead, AwxGridHead] = None
    palette: Optional[np.ndarray] = None
    calibration: Optional[np.ndarray] = None
    calibrating: Optional[bool] = False
    positioning: Optional[AwxPositioning] = None
    data: Union[np.ndarray, dict] = None
    missing_value: Optional[float] = None
    _ds: Optional[xr.DataArray] = None

    def __post_init__(self):
        if self.pathfile is not None:
            self.read(self.pathfile)

    def sel(self, lat=None, lon=None, method=None, **kwargs):
        """interface to select variable from file by given more filter and clip parameters

        Parameters
        ----------
        lat (slice, optional): latitude extent
        lon (slice, optional): longitude extent
        method (str): method in processing un-precise boundary, see `xarray.DataArray.sel`
        kwargs (dict): other parameters used to process netcdf variable

        Returns
        -------
        ncvar (xarray.DataArray): Readed variable in xarray.DataArray
        """
        darray = self.clip(self._ds, lat, lon, method)
        ncvar = darray.sel(**kwargs)

        return ncvar

    def clip(self, darray, lat=None, lon=None, method=None):
        """Clip the longitude and latiutde of AWX DataArray

        Parameters
        ----------
        darray (xarray.DataArray): DataArray with longitude and latitude coordinates
        lat (slice, optional): latitude extent
        lon (slice, optional): longitude extent
        method (str): method in processing un-precise boundary, see `xarray.DataArray.sel`

        Returns
        -------
        darray (xarray.DataArray): regional clipped DataArray
        """
        if lat is not None or lon is not None:
            if len(darray['lat'].shape) == 2:
                if isinstance(lat, slice):
                    start, stop = min(lat.start, lat.stop), max(lat.start, lat.stop)
                    darray = darray.where((darray['lat'] >= start) & (darray['lat'] <= stop),
                                          drop=True)
                if isinstance(lon, slice):
                    darray = darray.where(
                        (darray['lon'] >= lon.start) & (darray['lon'] <= lon.stop),
                        drop=True)
            else:
                if lat is not None:
                    if darray.indexes['lat'].is_monotonic_decreasing:
                        lat = slice(lat.stop, lat.start, lat.step)
                    darray = darray.sel({'lat': lat}, method=method)
                if lon is not None:
                    darray = darray.sel({'lon': lon}, method=method)
        return darray

    def read(self, path_or_bytes: Union[str, PathLike, bytes]):
        """read AWX file or bytes

        Parameters
        ----------
        path_or_bytes: path to AWX file or data in bytes
        """
        if isinstance(path_or_bytes, (str, PathLike)):
            self._deconstruct_filename(path_or_bytes)
            with open(self.pathfile, 'rb') as f:
                bytes_array = f.read()
        else:
            bytes_array = path_or_bytes

        head1_type = AwxHead.dtype()
        p = calcsize(head1_type)
        head1_data = unpack(head1_type, bytes_array[:p])
        self.head1 = AwxHead(*head1_data)

        Head2 = self.head2_interface()
        if Head2 is None:
            raise NotImplementedError("Unsupported file format.")
        head2_type = Head2.dtype()
        p_n = p + calcsize(head2_type)
        head2_data = unpack(head2_type, bytes_array[p:p_n])
        p = p_n
        self.head2 = Head2(*head2_data)

        pp = self.head1.head_record * self.head1.record_length
        if self.head1.product_kind in [1, 2]:
            if self.head2.palette_data_length > 0:
                p_n = p + self.head2.palette_data_length
                self.palette = np.frombuffer(bytes_array[p:p_n],
                                             dtype=np.uint8, count=self.head2.palette_data_length)
                p = p_n
            if self.head2.calibration_data_length > 0:
                p_n = p + self.head2.calibration_data_length
                self.calibration = np.frombuffer(bytes_array[p:p_n], dtype=np.uint16,
                                                 count=int(self.head2.calibration_data_length / 2))
                p = p_n
            if self.head2.position_data_length > 0:
                p_n = p + 16
                self.positioning = AwxPositioning(*unpack('8h', bytes_array[p:p_n]))
            n_obs = self.head2.width * self.head2.height
            data = np.frombuffer(bytes_array[pp: pp + n_obs], dtype=np.uint8, count=n_obs)
            if self.calibrating:
                data = self.calibration[data]
            data = np.array(data).reshape(self.head2.height, self.head2.width)
        elif self.head1.product_kind == 3:
            n_obs = self.head2.size_h * self.head2.size_v
            body_data = np.frombuffer(bytes_array[pp: pp + n_obs * self.head2.grid_byte],
                                      dtype=f'u{self.head2.grid_byte}', count=n_obs).astype('f')
            if self.head2.has_quality:
                down, up = self.head2.quality_down, self.head2.quality_up
                body_data[(body_data < down) | (body_data > up)] = np.nan
            data = (body_data + self.head2.grid_base) / self.head2.grid_scale
            data = np.array(data).reshape(self.head2.size_v, self.head2.size_h)
        else:
            raise NotImplementedError("Not supported product kind.")
        self.data = data
        self._ds = self.to_xarray()

    def to_xarray(self) -> xr.DataArray:
        """convert to xarray.DataArray"""
        attrs = self.head1.__dict__
        attrs.update(self.head2.__dict__)
        if self.head1.product_kind in (1, 2):
            coords = self.image_coordinate()
            if 'x' in coords and 'y' in coords:
                attrs['proj4'] = coords.pop('proj4')
                return xr.DataArray(self.data[None, :, :], dims=('time', 'y', 'x'),
                                    coords=coords, attrs=attrs)
            else:
                return xr.DataArray(self.data[None, :, :],
                                    dims=('time', 'lat', 'lon'),
                                    coords=coords, attrs=attrs)
        elif self.head1.product_kind == 3:
            attrs['description'] = self.head2.grid_elements_dict()[self.head2.grid_element]
            return xr.DataArray(self.data[None, :, :], dims=('time', 'lat', 'lon'),
                                coords=self.grid_coordinate(), attrs=attrs)
        else:
            raise NotImplementedError("Not supported product kind.")

    def grid_coordinate(self) -> dict:
        """calculate AWX grid field coordination"""
        h2 = self.head2
        time = datetime(h2.start_year, h2.start_month, h2.start_day, h2.start_hour, h2.start_minute)
        lats = np.linspace(h2.ul_lat / 100., h2.lr_lat / 100., h2.size_v)
        lons = np.linspace(h2.ul_lon / 100., h2.lr_lon / 100., h2.size_h)
        return {'time': [time], 'lat': lats, 'lon': lons}

    def image_coordinate(self) -> dict:
        """calculate AWX image field coordination"""
        h2 = self.head2
        if hasattr(h2, 'year'):
            time = datetime(h2.year, h2.month, h2.day, h2.hour, h2.minute)
        else:
            time = datetime(h2.start_year, h2.start_month, h2.start_day,
                            h2.start_hour, h2.start_minute)
        proj_str = self.proj_str()
        if proj_str is None:
            lats = np.linspace(h2.ul_lat / 100., h2.lr_lat / 100., h2.height)
            lons = np.linspace(h2.ul_lon / 100., h2.lr_lon / 100., h2.width)
            return {'time': [time], 'lat': lats, 'lon': lons}
        else:
            proj = CRS.from_proj4(proj_str)
            transformer = Transformer.from_crs("EPSG:4326", proj, always_xy=True)
            cx, cy = transformer.transform(h2.clon / 100., h2.clat / 100.)
            dx, dy = h2.reso_h / 100. * 1000, h2.reso_v / 100. * 1000,
            ll_x = cx - (dx * h2.width / 2.)
            ll_y = cy - (dy * h2.height / 2.)
            ur_x = cx + (dx * (h2.width / 2. - 1))
            ur_y = cy + (dy * (h2.height / 2. - 1))
            x = np.linspace(ll_x, ur_x, h2.width)
            y = np.linspace(ur_y, ll_y, h2.height)

            transformer = Transformer.from_crs(proj, "EPSG:4326", always_xy=True)
            x2d, y2d = np.meshgrid(x, y)
            lons, lats = transformer.transform(x2d, y2d)
            lons = xr.DataArray(lons, dims=('y', 'x'), coords={'y': y, 'x': x})
            lats = xr.DataArray(lats, dims=('y', 'x'), coords={'y': y, 'x': x})
            return {'time': [time], 'lat': lats, 'lon': lons, 'x': x, 'y': y, 'proj4': proj_str}

    def proj_scale_factor(self, std_parallel):
        """calculate projection scale factor"""
        from math import pi, sin, cos, tan, pow

        e = 0.081819191
        std_parallel = std_parallel * pi / 180.
        if std_parallel > 0:
            tf = tan(pi / 4. - std_parallel / 2.) * (
                pow((1. + e * sin(std_parallel)) / (1. - e * sin(std_parallel)), e / 2.))
        else:
            tf = tan(pi / 4. + std_parallel / 2.) * (
                pow((1. + e * sin(std_parallel)) / (1. - e * sin(std_parallel)), e / 2.))
        mf = cos(std_parallel) / pow(1. - e * e * pow(sin(std_parallel), 2.), 0.5)
        k0 = mf * (pow(pow(1. + e, 1. + e) * pow(1. - e, 1. - e), 0.5)) / (2. * tf)
        return k0

    def proj_str(self):
        """generate projection string to construct Proj"""
        h2 = self.head2
        proj_type = h2.projection
        if proj_type == 1:
            proj = f'+proj=lcc +lon_0={h2.clon / 100.} +lat_0={h2.clat / 100.} ' \
                   f'+lat_1={h2.std_lat1_or_lon / 100.} +lat_2={h2.std_lat2 / 100.}'
        elif proj_type == 2:
            # TODO: find the reason why the std lat in data head is wrong ?
            # set std lat to 0 (default), otherwise the reprojected coordination error
            # proj = f'+proj=merc +lon_0={h2.clon / 100.} +lat_ts={h2.std_lat1_or_lon / 100.}'
            proj = f'+proj=merc +lon_0={h2.clon / 100.}'
        elif proj_type == 3:
            proj = f'+proj=stere +lat_0={h2.clat / 100.} +lon_0={h2.clon / 100.} ' \
                   f'+k={self.proj_scale_factor(h2.std_lat1_or_lon / 100.)}'
        elif proj_type == 4:
            proj = None
        else:
            raise NotImplementedError("Not supported projection type.")
        return proj

    def head2_interface(self):
        kind = self.head1.product_kind
        heads = [None, AwxGeosImageHead, AwxPolarImageHead, AwxGridHead, None, None]
        return heads[kind]

    def _deconstruct_filename(self, filename):
        """deconstruct AWX filename to get base information of file"""
        name_segs = os.path.splitext(os.path.basename(filename))[0].split('_')
        if len(name_segs) == 6:
            self.name = AwxFileName(*name_segs)
        elif len(name_segs) == 12:
            self.name = AwxDaasFileName(*name_segs)
        else:
            pass

    @property
    def values(self) -> xr.DataArray:
        return self._ds

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close open dataset"""
        del self.data, self.head1, self.head2, self.palette, self.calibration
        del self.positioning, self._ds, self.pathfile

    def __repr__(self) -> str:
        """print"""
        return self.head1.__repr__() + '\n' + self.head2.__repr__()


def _info():
    example_text = """Example:
     awx_info FY2G_TBB_IR1_OTG_20150729_0000.AWX
     """
    parser = argparse.ArgumentParser(description='Model Data Dumper', epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='awx file name')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    awx = Awx(args.filename)
    print(awx.head1, '\n', awx.head2, '\n', awx._ds)


def _convert_to_nc():
    example_text = """Example:
     awx_to_nc FY2G_TBB_IR1_OTG_20150729_0000.AWX FY2G_TBB_IR1_OTG_20150729_0000.nc
     """
    parser = argparse.ArgumentParser(description='Model Data Dumper', epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='awx file name')
    parser.add_argument('outfile', help='path to convert netcdf file name')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    awx = Awx(args.filename)
    ds = awx.values
    ds.name = 'AWX'
    ds = ds.to_dataset()
    comp = dict(zlib=True, complevel=5)
    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(args.outfile, encoding=encoding)
