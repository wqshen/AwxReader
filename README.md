# 读取卫星中心9210格点场数据文件（后缀AWX）
本工程提供AWX文件的读取和格式转换以及快速绘图的功能。

## 例：

```python
    from awx import AwxGridField
    pathfile = r'../data/FY2E_CTA_MLT_OTG_20170126_0130.AWX'
    agf = AwxGridField(pathfile)
    # 输出文件信息
    print(agf)
    # 使用经纬度裁剪数据
    agf.clipper((slice(30,27), slice(118,121)))
    # 绘图
    agf.plot()
    # 输出netCDF4压缩格式（依赖于netCDF4库）
    # agf.to_netcdf('test.nc')
    # 转换为xarray.DataArray对象（依赖于xarray库）
    # da = agf.to_dataarray()
    # print(da)
    # 转换为pandas.DataFrame对象（依赖于Pandas库）
    # df = agf.to_dataframe()
    # print(df)
