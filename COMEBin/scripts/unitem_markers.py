# Adapted from https://github.com/dparks1134/UniteM/blob/master/unitem/markers.py
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import logging

from scripts.unitem_defaults import *


class Markers():
    """Read bin quality from CheckM2 output."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def read_quality_report(self, quality_report_tsv):
        """Read per-bin quality from a CheckM2 quality_report.tsv file.

        Parameters
        ----------
        quality_report_tsv : str
            Path to CheckM2's quality_report.tsv.

        Returns
        -------
        dict
            Mapping of bin name (str) to (completeness, contamination) tuple.
        """
        quality = {}
        if not os.path.exists(quality_report_tsv):
            self.logger.warning('CheckM2 quality report not found: %s' % quality_report_tsv)
            return quality

        with open(quality_report_tsv) as f:
            f.readline()  # skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    name = parts[0]
                    try:
                        comp = float(parts[1])
                        cont = float(parts[2])
                        quality[name] = (comp, cont)
                    except ValueError:
                        continue
        return quality
