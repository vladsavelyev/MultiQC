""" MultiQC module to parse output from QUAST """


import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="My Module",
            anchor="mymodule",
            href="http://mymodule.ru/",
            info="is my awesome module.",
            doi="10.1093/bioinformatics/btt086",
        )

        self.data = self.ignore_samples(self.data)

        self.data = {}
        if len(self.data) == 0:
            raise UserWarning

        self.write_data_file(self.data, "multiqc_mymodule")