"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.meas.extensions.astrometryNet


_g = globals()
_g.update(build_package_configs(
    project_name='meas_extensions_astrometryNet',
    version=lsst.meas.extensions.astrometryNet.version.__version__))
