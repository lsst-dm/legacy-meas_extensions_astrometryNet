from __future__ import absolute_import, division, print_function

from lsst.utils import continueClass

from . import MultiIndex


@continueClass  # noqa
class MultiIndex:
    def __iter__(self):
        """Get an iterator over the indices in a MultiIndex

        Do not modify the number or location of the indices while using the iterator.
        """
        for i in range(len(self)):
            yield self[i]
