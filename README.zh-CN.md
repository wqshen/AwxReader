# 读取气象卫星分发产品（后缀AWX）

本工程提供AWX文件的读取，目前支持如下类型

- 图像产品
- 格点场定量产品

## 包的安装

使用 `pip` 安装

```shell
pip install awx
```


## 快速入门使用

### 包的基本使用方法

**1 读取文件、访问数据、经纬度切片、转存netCDF文件**

```python

from awx import Awx

pathfile = r'../data/ANI_VIS_R02_20230217_1000_FY2G.AWX'
ds = Awx(pathfile)

# 打印文件头信息
print(ds)

# 获取xarray.DataArray格式的观测数据
print(ds.values)

# 裁剪到指定的经纬度范围
print(ds.sel(lat=slice(20, 40), lon=slice(100, 130)))

# 保存数据为netCDF格式文件
ds.values.to_netcdf('ANI_VIS_R02_20230217_1000_FY2G.nc')
```

**2 利用Matplotlib绘制无投影的云图**

```python
import matplotlib.pyplot as plt
from awx import Awx

fpath = r'./data/ANI_VIS_R02_20230217_1000_FY2G.AWX'
ds = Awx(pathfile=fpath)
print(ds)
dar = ds.values.squeeze()
plt.pcolormesh(dar.lon, dar.lat, dar, cmap='Greys_r')
plt.savefig('ANI_VIS_R02_20230217_1000_FY2G_NoProj.png', dpi=300)
plt.show()
```

![ANI_VIS_R02_20230217_1000_FY2G_NoProj.png](doc%2FANI_VIS_R02_20230217_1000_FY2G_NoProj.png)

**3 以数据指定的投影方式绘制云图**

```python
import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from awx import Awx

# fpath = r'./data/ANI_VIS_R02_20230217_1000_FY2G.AWX'  # Mercator
fpath = r'./data/ANI_IR2_R01_20230217_0800_FY2G.AWX'  # lambert
ds = Awx(pathfile=fpath)
print(ds)
dar = ds.values.squeeze()

plt.figure(figsize=(8, 8))

if dar.projection == 1:
    proj = ccrs.LambertConformal(central_longitude=dar.clon / 100,
                                 central_latitude=dar.clat / 100,
                                 standard_parallels=(dar.std_lat1_or_lon / 100.,
                                                     dar.std_lat2 / 100.))
    extent = [dar.x.min(), dar.x.max(), dar.y.min(), dar.y.max()]
elif dar.projection == 2:
    proj = ccrs.Mercator(central_longitude=dar.clon / 100,
                         latitude_true_scale=dar.std_lat1_or_lon / 100.)
    extent = [dar.x.min(), dar.x.max(), dar.y.min(), dar.y.max()]
elif dar.projection == 4:
    proj = ccrs.PlateCarree(central_longitude=dar.clon / 100.)
    extent = [dar.lon.min(), dar.lon.max(), dar.lat.min(), dar.lat.max()]
else:
    raise NotImplementedError()
ax = plt.axes(projection=proj)
ax.set_extent(extent, crs=proj)
ax.coastlines(resolution='110m')
ax.gridlines(draw_labels=True)
ax.pcolormesh(dar.x, dar.y, dar, cmap='Greys_r')
plt.savefig(os.path.splitext(os.path.basename(fpath))[0] + '.png', dpi=300, bbox_inches='tight')
plt.show()

```

![ANI_VIS_R02_20230217_1000_FY2G.png](doc%2FANI_VIS_R02_20230217_1000_FY2G.png)

![ANI_IR2_R01_20230217_0800_FY2G.png](doc%2FANI_IR2_R01_20230217_0800_FY2G.png)

### 命令行程序

#### awx_info

打印AWX文件头信息

用法:
    
    awx_info AWX文件名

示例:

    `awx_info FY2G_TBB_IR1_OTG_20150729_0000.AWX`

#### awx_to_nc

转换AWX文件为netCDF格式

用法:

    awx_to_nc AWX文件名 NetCDF文件名

示例:

    `awx_to_nc FY2G_TBB_IR1_OTG_20150729_0000.AWX FY2G_TBB_IR1_OTG_20150729_0000.nc`