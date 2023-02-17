# -*- coding: utf-8 -*-
# @Author: wqshen
# @Email: wqshen91@gmail.com
# @Date: 2023/2/9 22:07
# @Last Modified by: wqshen

from typing import cast
from dataclasses import dataclass, fields


@dataclass
class AwxGridHead:
    satellite_name: str
    grid_element: int
    grid_byte: int
    grid_base: int
    grid_scale: int
    time_scale: int
    start_year: int
    start_month: int
    start_day: int
    start_hour: int
    start_minute: int
    end_year: int
    end_month: int
    end_day: int
    end_hour: int
    end_minute: int
    ul_lat: int
    ul_lon: int
    lr_lat: int
    lr_lon: int
    grid_unit: int
    reso_h: int
    reso_v: int
    size_h: int
    size_v: int
    has_land: int = 0
    land: int = 0
    has_cloud: int = 0
    cloud: int = 0
    has_water: int = 0
    water: int = 0
    has_ice: int = 0
    ice: int = 0
    has_quality: int = 0
    quality_up: int = 0
    quality_down: int = 0
    reserve: int = 0

    @staticmethod
    def dtype() -> str:
        return '=8s36h'

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            if isinstance(field_value, bytes):
                field_value = field_value.rstrip(b'\x00').decode('utf8')
            setattr(self, field.name, field.type(field_value))

    def grid_elements_dict(self):
        return {0: 'numerical weather prediction',
                1: 'sea surface temperature (K) ',
                2: 'sea ice distribution (dimensionless) ',
                3: 'sea ice density (dimensionless)',
                4: 'outgoing longwave radiation (W/m2)',
                5: 'normalized vegetation index (dimensionless)',
                6: 'vegetation index ratio (dimensionless) ',
                7: 'snow distribution (dimensionless)',
                8: 'soil moisture (kg/m3)',
                9: 'sunshine (hour) ',
                10: 'cloud top height (hPa) ',
                11: 'cloud top temperature (K)',
                12: 'low cirrus (dimensionless)',
                13: 'high cirrus (dimensionless)',
                14: 'precipitation index (mm/1h) ',
                15: 'precipitation index (mm/6h)',
                16: 'precipitation index (mm/12h) ',
                17: 'precipitation index (mm/24h) ',
                18: 'upper tropospheric water vapor (relative humidity) (dimensionless)',
                19: 'brightness temperature',
                20: 'cloud amount (percentage)',
                21: 'cloud classification (dimensionless)',
                22: 'precipitation estimates (mm/6h)',
                23: 'precipitation estimates (mm/24h)',
                24: 'clear sky atmospheric precipitation (mm)',
                25: 'reserved',
                26: 'ground incident solar radiation (W/m2)',
                27: 'reserved ',
                28: 'reserved ',
                29: 'reserved ',
                30: 'reserved ',
                31: 'cloud humidity profiles (1000 hPa)',
                32: 'cloud humidity profiles (850 hPa)',
                33: 'cloud humidity profiles (700 hPa)',
                34: 'cloud humidity profiles (600 hPa)',
                35: 'cloud humidity profiles (500 hPa)',
                36: 'cloud humidity profiles (400 hPa)',
                37: 'cloud humidity profiles (300 hPa)',
                101: 'clear sky environment monitoring datasets',
                201: 'ATOVS (1000-10 hPa) temperature fields (K)',
                202: 'ATOVS (1000-10 hPa) temperature fields (K)',
                203: 'ATOVS (1000-10 hPa) temperature fields (K)',
                204: 'ATOVS (1000-10 hPa) temperature fields (K)',
                205: 'ATOVS (1000-10 hPa) temperature fields (K)',
                206: 'ATOVS (1000-10 hPa) temperature fields (K)',
                207: 'ATOVS (1000-10 hPa) temperature fields (K)',
                208: 'ATOVS (1000-10 hPa) temperature fields (K)',
                209: 'ATOVS (1000-10 hPa) temperature fields (K)',
                210: 'ATOVS (1000-10 hPa) temperature fields (K)',
                211: 'ATOVS (1000-10 hPa) temperature fields (K)',
                212: 'ATOVS (1000-10 hPa) temperature fields (K)',
                213: 'ATOVS (1000-10 hPa) temperature fields (K)',
                214: 'ATOVS (1000-10 hPa) temperature fields (K)',
                215: 'ATOVS (1000-10 hPa) temperature fields (K)',

                401: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',
                402: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',
                403: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',
                404: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',
                405: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',
                406: 'ATOVS (1000 ~ 300 hPa) dew point temperature fields (K)',

                501: 'ATOVS atmospheric stability index (dimensionless)',
                502: 'ATOVS clear sky total atmospheric column water vapor content (mm)',
                503: 'ATOVS total atmospheric column ozone content (Db)',
                504: ' ATOVS outgoing longwave radiation (W/m2)',
                505: 'ATOVS cloud top height (hPa) ',
                506: 'ATOVS cloud top temperature (K)',
                507: 'ATOVS cloudiness (dimensionless)(ZK)',
                }
