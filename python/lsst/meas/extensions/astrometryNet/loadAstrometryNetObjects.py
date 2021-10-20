__all__ = ["LoadAstrometryNetObjectsTask", "LoadAstrometryNetObjectsConfig"]

import lsst.pipe.base as pipeBase
from lsst.meas.algorithms import LoadReferenceObjectsTask, getRefFluxField
from lsst.meas.algorithms.loadReferenceObjects import convertToNanojansky
from lsst.utils.timer import timeMethod
from . import astrometry_net
from .multiindex import AstrometryNetCatalog, getConfigFromEnvironment

LoadAstrometryNetObjectsConfig = LoadReferenceObjectsTask.ConfigClass

# The following block adds links to this task from the Task Documentation page.
## \addtogroup LSST_task_documentation
## \{
## \page measAstrom_loadAstrometryNetObjectsTask
## \ref LoadAstrometryNetObjectsTask "LoadAstrometryNetObjectsTask"
##      Load reference objects from astrometry.net index files
## \}


class LoadAstrometryNetObjectsTask(LoadReferenceObjectsTask):
    """!Load reference objects from astrometry.net index files

    @anchor LoadAstrometryNetObjectsTask_

    @section meas_astrom_loadAstrometryNetObjects_Contents Contents

     - @ref meas_astrom_loadAstrometryNetObjects_Purpose
     - @ref meas_astrom_loadAstrometryNetObjects_Initialize
     - @ref meas_astrom_loadAstrometryNetObjects_IO
     - @ref meas_algorithms_loadReferenceObjects_Schema
     - @ref meas_astrom_loadAstrometryNetObjects_Config
     - @ref meas_astrom_loadAstrometryNetObjects_Example
     - @ref meas_astrom_loadAstrometryNetObjects_Debug

    @section meas_astrom_loadAstrometryNetObjects_Purpose  Description

    Load reference objects from astrometry.net index files.

    @section meas_astrom_loadAstrometryNetObjects_Initialize   Task initialisation

    @copydoc \_\_init\_\_

    @section meas_astrom_loadAstrometryNetObjects_IO       Invoking the Task

    @copydoc loadObjectsInBBox

    @section meas_astrom_loadAstrometryNetObjects_Config       Configuration parameters

    See @ref LoadAstrometryNetObjectsConfig

    @section meas_astrom_loadAstrometryNetObjects_Example  A complete example of using
        LoadAstrometryNetObjectsTask

    LoadAstrometryNetObjectsTask is a subtask of AstrometryTask, which is called by PhotoCalTask.
    See \ref pipe_tasks_photocal_Example.

    @section meas_astrom_loadAstrometryNetObjects_Debug        Debug variables

    LoadAstrometryNetObjectsTask does not support any debug variables.
    """
    ConfigClass = LoadAstrometryNetObjectsConfig

    def __init__(self, config=None, andConfig=None, **kwargs):
        """!Create a LoadAstrometryNetObjectsTask

        @param[in] config  configuration (an instance of self.ConfigClass); if None use self.ConfigClass()
        @param[in] andConfig  astrometry.net data config (an instance of AstromNetDataConfig, or None);
            if None then use andConfig.py in the astrometry_net_data product (which must be setup)
        @param[in] kwargs  additional keyword arguments for pipe_base Task.\_\_init\_\_

        @throw RuntimeError if andConfig is None and the configuration cannot be found,
            either because astrometry_net_data is not setup in eups
            or because the setup version does not include the file "andConfig.py"
        """
        LoadReferenceObjectsTask.__init__(self, config=config, **kwargs)
        self.andConfig = andConfig
        self.haveIndexFiles = False  # defer reading index files until we know they are needed
        # because astrometry may not be used, in which case it may not be properly configured

    @timeMethod
    def loadSkyCircle(self, ctrCoord, radius, filterName=None, epoch=None, centroids=True):
        """!Load reference objects that overlap a circular sky region

        @param[in] ctrCoord  center of search region (an afwGeom.Coord)
        @param[in] radius  radius of search region (an geom.Angle)
        @param[in] filterName  name of filter, or None for the default filter;
            used for flux values in case we have flux limits (which are not yet implemented)
        @param[in] epoch  Epoch for proper motion and parallax correction
                    (an astropy.time.Time), or None
        centroids : `bool` (optional)
            Ignored: a.net refcats always have centroid fields.

        No proper motion correction is made, since our astrometry.net catalogs
        typically don't support that, and even if they do they format is uncertain.
        Users interested in proper motion corrections should use the
        lsst.meas.algorithms.LoadIndexedReferenceObjectsTask or they will need to
        subclass and define how the proper motion correction is to be done.

        @return an lsst.pipe.base.Struct containing:
        - refCat a catalog of reference objects with the
            \link meas_algorithms_loadReferenceObjects_Schema standard schema \endlink
            as documented in LoadReferenceObjects, including photometric, resolved and variable;
            hasCentroid is False for all objects.
        - fluxField = name of flux field for specified filterName
        """
        self._readIndexFiles()

        names = []
        mcols = []
        ecols = []
        for col, mcol in self.andConfig.magColumnMap.items():
            names.append(col)
            mcols.append(mcol)
            ecols.append(self.andConfig.magErrorColumnMap.get(col, ''))
        margs = (names, mcols, ecols)

        solver = self._getSolver()

        # Find multi-index files within range
        multiInds = self._getMIndexesWithinRange(ctrCoord, radius)

        # compute solver.getCatalog arguments that follow the list of star kd-trees:
        # - center equatorial angle (e.g. RA) in deg
        # - center polar angle (e.g. Dec) in deg
        # - radius, in deg
        # - idColumn
        # - (margs)
        # - star-galaxy column
        # - variability column
        fixedArgTuple = (
            ctrCoord,
            radius,
            self.andConfig.idColumn,
        ) + margs + (
            self.andConfig.starGalaxyColumn,
            self.andConfig.variableColumn,
            True,  # eliminate duplicate IDs
        )

        self.log.debug("search for objects at %s with radius %s deg", ctrCoord, radius.asDegrees())
        with LoadMultiIndexes(multiInds):
            # We just want to pass the star kd-trees, so just pass the
            # first element of each multi-index.
            inds = tuple(mi[0] for mi in multiInds)
            refCat = solver.getCatalog(inds, *fixedArgTuple)

        self._addFluxAliases(schema=refCat.schema)

        fluxField = getRefFluxField(schema=refCat.schema, filterName=filterName)

        # NOTE: sourceSelectors require contiguous catalogs, so ensure
        # contiguity now, so views are preserved from here on.
        if not refCat.isContiguous():
            refCat = refCat.copy(deep=True)

        # Update flux fields to be nJy. a.net catalogs do not have a conversion script.
        self.log.warn("Loading A.net reference catalog with old style units in schema.")
        self.log.warn("A.net reference catalogs will not be supported in the future.")
        self.log.warn("See RFC-562 and RFC-575 for more details.")
        refCat = convertToNanojansky(refCat, self.log)

        self.log.debug("found %d objects", len(refCat))
        return pipeBase.Struct(
            refCat=refCat,
            fluxField=fluxField,
        )

    @timeMethod
    def _readIndexFiles(self):
        """!Read all astrometry.net index files, if not already read
        """
        if self.haveIndexFiles:
            return

        self.log.debug("read index files")
        self.haveIndexFiles = True  # just try once

        if self.andConfig is None:
            self.andConfig = getConfigFromEnvironment()

        self.multiInds = AstrometryNetCatalog(self.andConfig)

    def _getMIndexesWithinRange(self, ctrCoord, radius):
        """!Get list of muti-index objects within range

        @param[in] ctrCoord  center of search region (an afwGeom.Coord)
        @param[in] radius  radius of search region (an geom.Angle)

        @return list of multiindex objects
        """
        return [mi for mi in self.multiInds if mi.isWithinRange(ctrCoord, radius)]

    def _getSolver(self):
        solver = astrometry_net.Solver()
        # HACK, set huge default pixel scale range.
        lo, hi = 0.01, 3600.
        solver.setPixelScaleRange(lo, hi)
        return solver


class LoadMultiIndexes:
    """Context manager for loading and unloading astrometry.net multi-index files
    """

    def __init__(self, multiInds):
        self.multiInds = multiInds

    def __enter__(self):
        for mi in self.multiInds:
            mi.reload()
        return self.multiInds

    def __exit__(self, typ, val, trace):
        for mi in self.multiInds:
            mi.unload()
