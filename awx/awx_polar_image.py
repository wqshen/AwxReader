# -*- coding: utf-8 -*-
# @Author: wqshen
# @Email: wqshen91@gmail.com
# @Date: 2023/2/9 22:07
# @Last Modified by: wqshen

from typing import cast
from dataclasses import dataclass, fields


@dataclass
class AwxPolarImageHead:
    satellite_name: str
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
    channel: int
    channel_r: int
    channel_g: int
    channel_b: int
    track_mark: int
    track_id: int
    pixel_bytes: int
    projection: int
    product_type: int
    width: int
    height: int
    ul_row_id: int
    ul_pixel_id: int
    sampling_rate: int
    ul_lat: int
    lr_lat: int
    ul_lon: int
    lr_lon: int
    clat: int
    clon: int
    std_lat1_or_lon: int
    std_lat2: int
    reso_h: int
    reso_v: int
    mark_geogrid_overlap: int
    value_geogrid_overlap: int
    palette_data_length: int = 0
    calibration_data_length: int = 0
    position_data_length: int = 0
    reserved: int = 0

    @staticmethod
    def dtype() -> str:
        return '=8s40h'

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            if isinstance(field_value, bytes):
                field_value = field_value.rstrip(b'\x00').decode('utf8')
            setattr(self, field.name, field.type(field_value))
