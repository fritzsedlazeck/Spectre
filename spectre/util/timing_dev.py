from spectre.util import logger
from datetime import datetime
import numpy as np

class DevTiming(object):
    def __init__(self, timing_process=None, name=None):
        self.start_ts = None
        self.timing_process = timing_process if timing_process is not None else "process"
        self.logger = logger.setup_log(name, True) if name is not None else logger.setup_log(__name__, True)

    def start(self):
        # Output
        self.logger.debug(f'Starting {self.timing_process} ...')
        # time start for process
        start_dt = datetime.now()
        self.start_ts = datetime.timestamp(start_dt)

    def end(self, add_sep=True, return_time=False):
        end_dt = datetime.now()
        end_ts = datetime.timestamp(end_dt)
        separator = "# ########################## #"
        time_spent = np.round(end_ts - self.start_ts, 3)
        self.logger.debug(f'Timing: {time_spent}\n{separator}') if add_sep \
            else self.logger.debug(f'Timing: {time_spent}')
        if return_time:
            return time_spent
