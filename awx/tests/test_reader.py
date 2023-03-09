# -*- coding: utf-8 -*-
# @Author: wqshen
# @Email: wqshen91@gmail.com
# @Date: 2021/5/31 15:47
# @Last Modified by: wqshen

import os
import pytest
from awx import Awx


class TestReader:

    def test_read_product_3_fy2g_tbb(self):
        fpath = r'../../sampledata/FY2G_TBB_IR1_OTG_20150729_0000.AWX'
        ds = Awx(pathfile=fpath)
        print(ds)
        print(ds.values)
        print(ds.sel(lat=slice(20, 40), lon=slice(100, 130)))

    def test_read_product_3_fy2e_cta(self):
        fpath = r'../../sampledata/FY2E_CTA_MLT_OTG_20170126_0130.AWX'
        ds = Awx(pathfile=fpath)
        print(ds)
        print(ds.values)
        print(ds.sel(lat=slice(20, 40), lon=slice(100, 130)))

    def test_read_product_1_fy2_ir(self):
        fpath = r'../../sampledata/ANI_IR2_R01_20230217_0800_FY2G.AWX'
        ds = Awx(pathfile=fpath)
        print(ds)
        print(ds.values)
        print(ds.sel(lat=slice(20, 40), lon=slice(100, 130)))

    def test_read_product_1_fy2_vis(self):
        fpath = r'../../sampledata/ANI_VIS_R02_20230217_1000_FY2G.AWX'
        ds = Awx(pathfile=fpath)
        print(ds)
        print(ds.values)
        print(ds.sel(lat=slice(20, 40), lon=slice(100, 130)))

    def test_plot_matplotlib(self):
        import matplotlib.pyplot as plt

        fpath = r'../../sampledata/ANI_VIS_R02_20230308_1400_FY2G.AWX'
        ds = Awx(pathfile=fpath)
        print(ds)
        dar = ds.values.squeeze()
        plt.pcolormesh(dar.lon, dar.lat, dar, cmap='Greys_r')
        plt.savefig('ANI_VIS_R02_20230308_1400_FY2G_NoProj.png', dpi=300)
        plt.show()

    def get_proj_and_extent(self, dar):
        import cartopy.crs as ccrs

        if dar.projection == 1:
            proj = ccrs.LambertConformal(central_longitude=dar.clon / 100,
                                         central_latitude=dar.clat / 100,
                                         standard_parallels=(
                                             dar.std_lat1_or_lon / 100.,
                                             dar.std_lat2 / 100.
                                         ))
            extent = [dar.x.min(), dar.x.max(), dar.y.min(), dar.y.max()]
        elif dar.projection == 2:
            proj = ccrs.Mercator(central_longitude=dar.clon / 100,
                                 # latitude_true_scale=dar.std_lat1_or_lon / 100.
                                 )

            extent = [dar.x.min(), dar.x.max(), dar.y.min(), dar.y.max()]
        elif dar.projection == 4:
            proj = ccrs.PlateCarree()
            extent = [dar.lon.min(), dar.lon.max(), dar.lat.min(), dar.lat.max()]
        else:
            raise NotImplementedError()
        return proj, extent

    def test_plot_cartopy(self):
        import matplotlib.pyplot as plt

        # fpath = r'E:\ANI_VIS_R01_20230308_1400_FY2G.AWX'  # lambert VIS
        fpath = r'E:\ANI_VIS_R02_20230308_1400_FY2G.AWX'  # Mercator VIS
        ds = Awx(pathfile=fpath)
        print(ds)
        dar = ds.values.squeeze()
        proj, extent = self.get_proj_and_extent(dar)

        plt.figure(figsize=(8, 8))
        ax = plt.axes(projection=proj)
        ax.set_extent(extent, crs=proj)
        ax.coastlines(resolution='50m')
        ax.gridlines(draw_labels=True)
        if dar.projection == 4:
            ax.pcolormesh(dar.lon, dar.lat, dar, cmap='Greys_r')
        else:
            ax.pcolormesh(dar.x, dar.y, dar, cmap='Greys_r')
        plt.savefig(os.path.splitext(os.path.basename(fpath))[0] + '.png', dpi=300,
                    bbox_inches='tight')
        plt.show()


if __name__ == '__main__':
    pytest.main(['-q', 'test_diamond_reader.py'])
