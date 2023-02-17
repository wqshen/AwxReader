# -*- coding: utf-8 -*-
# @Author: wqshen
# @Email: wqshen91@gmail.com
# @Date: 2023/2/9 22:07
# @Last Modified by: wqshen

# Note: Not finished.


from typing import cast
from dataclasses import dataclass, fields


@dataclass
class AwxDiscreteHead:
    satellite_name: str
    element: int
    chars_by_record: int
    points: int
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
    inversion_method: int
    initial_field: int
    default: int = 0

    @staticmethod
    def dtype() -> str:
        return '=8s16h'

    def __post_init__(self):
        self._cast_fields_types()

    def _cast_fields_types(self):
        for field in fields(cast(dataclass, self)):
            field_value = getattr(self, field.name)
            if isinstance(field_value, bytes):
                field_value = field_value.rstrip(b'\x00').decode('utf8')
            setattr(self, field.name, field.type(field_value))


@dataclass
class AwxAtovs:
    pass


@dataclass
class AwxWind:
    pass
