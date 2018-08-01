#
# LSST Data Management System
# Copyright 2008-2017 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os.path
import math
import unittest

import lsst.utils.tests
import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.meas.base as measBase
from lsst.meas.extensions.astrometryNet import AstrometryNetDataConfig, \
    ANetAstrometryTask, LoadAstrometryNetObjectsTask
from test_findAstrometryNetDataDir import setupAstrometryNetDataDir


class TestAstrometricSolver(lsst.utils.tests.TestCase):

    def setUp(self):
        self.datapath = setupAstrometryNetDataDir('photocal')

        self.bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), lsst.geom.Extent2I(3001, 3001))
        crpix = lsst.geom.Box2D(self.bbox).getCenter()
        self.tanWcs = afwGeom.makeSkyWcs(crpix=crpix,
                                         crval=lsst.geom.SpherePoint(215.5, 53.0, lsst.geom.degrees),
                                         cdMatrix=afwGeom.makeCdMatrix(scale=5.1e-5*lsst.geom.degrees))
        self.exposure = afwImage.ExposureF(self.bbox)
        self.exposure.setWcs(self.tanWcs)
        self.exposure.setFilter(afwImage.Filter("r", True))
        andConfig = AstrometryNetDataConfig()
        andConfig.load(os.path.join(self.datapath, 'andConfig2.py'))
        andConfig.magErrorColumnMap = {}
        self.refObjLoader = LoadAstrometryNetObjectsTask(andConfig=andConfig)

    def tearDown(self):
        del self.tanWcs
        del self.exposure
        del self.refObjLoader

    def testTrivial(self):
        """Test fit with no distortion
        """
        self.doTest(afwGeom.makeIdentityTransform())

    def testRadial(self):
        """Test fit with radial distortion

        The offset comes from the fact that the CCD is not centered
        """
        self.doTest(afwGeom.makeRadialTransform([0, 1.01, 1e-7]))

    def makeSourceSchema(self):
        schema = afwTable.SourceTable.makeMinimalSchema()
        measBase.SingleFrameMeasurementTask(schema=schema)  # expand the schema
        return schema

    def doTest(self, pixelsToTanPixels):
        """Test using pixelsToTanPixels to distort the source positions
        """
        distortedWcs = afwGeom.makeModifiedWcs(pixelTransform=pixelsToTanPixels, wcs=self.tanWcs,
                                               modifyActualPixels=False)
        self.exposure.setWcs(distortedWcs)
        sourceCat = self.makeSourceCat(distortedWcs)
        print("number of stars =", len(sourceCat))
        config = ANetAstrometryTask.ConfigClass()
        solver = ANetAstrometryTask(config=config, refObjLoader=self.refObjLoader,
                                    schema=self.makeSourceSchema())
        results = solver.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        fitWcs = self.exposure.getWcs()
        self.assertRaises(Exception, self.assertWcsAlmostEqualOverBBox, fitWcs, distortedWcs)
        self.assertWcsAlmostEqualOverBBox(distortedWcs, fitWcs, self.bbox,
                                          maxDiffSky=0.01*lsst.geom.arcseconds, maxDiffPix=0.02)

        srcCoordKey = afwTable.CoordKey(sourceCat.schema["coord"])
        refCoordKey = afwTable.CoordKey(results.refCat.schema["coord"])
        maxAngSep = 0*lsst.geom.radians
        maxPixSep = 0
        for refObj, src, d in results.matches:
            refCoord = refObj.get(refCoordKey)
            refPixPos = fitWcs.skyToPixel(refCoord)
            srcCoord = src.get(srcCoordKey)
            srcPixPos = src.getCentroid()

            angSep = refCoord.separation(srcCoord)
            maxAngSep = max(maxAngSep, angSep)

            pixSep = math.hypot(*(srcPixPos-refPixPos))
            maxPixSep = max(maxPixSep, pixSep)
        print("max angular separation = %0.4f arcsec" % (maxAngSep.asArcseconds(),))
        print("max pixel separation = %0.3f" % (maxPixSep,))
        self.assertLess(maxAngSep.asArcseconds(), 0.005)
        self.assertLess(maxPixSep, 0.03)

        # try again, but without fitting the WCS
        config.forceKnownWcs = True
        solverNoFit = ANetAstrometryTask(config=config, refObjLoader=self.refObjLoader,
                                         schema=self.makeSourceSchema())
        self.exposure.setWcs(distortedWcs)
        noFitResults = solverNoFit.run(
            sourceCat=sourceCat,
            exposure=self.exposure,
        )
        self.assertGreater(len(noFitResults.refCat), 300)

    def makeSourceCat(self, distortedWcs):
        """Make a source catalog by reading the position reference stars and distorting the positions
        """
        loadRes = self.refObjLoader.loadPixelBox(bbox=self.bbox, wcs=distortedWcs, filterName="r")
        refCat = loadRes.refCat
        refCentroidKey = afwTable.Point2DKey(refCat.schema["centroid"])
        refFluxRKey = refCat.schema["r_flux"].asKey()

        sourceSchema = self.makeSourceSchema()
        sourceCat = afwTable.SourceCatalog(sourceSchema)
        sourceCentroidKey = afwTable.Point2DKey(sourceSchema["slot_Centroid"])
        sourceFluxKey = sourceSchema["slot_PsfFlux_flux"].asKey()
        sourceFluxErrKey = sourceSchema["slot_PsfFlux_fluxErr"].asKey()

        sourceCat.reserve(len(refCat))
        for refObj in refCat:
            src = sourceCat.addNew()
            src.set(sourceCentroidKey, refObj.get(refCentroidKey))
            src.set(sourceFluxKey, refObj.get(refFluxRKey))
            src.set(sourceFluxErrKey, refObj.get(refFluxRKey)/100)
        return sourceCat


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
