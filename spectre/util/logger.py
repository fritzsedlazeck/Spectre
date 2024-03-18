import sys
import logging

def setup_log(name, debug=False):
    py_logger = logging.getLogger(name)
    log_level = logging.INFO if not debug else logging.DEBUG
    py_logger.setLevel(log_level)
    py_formatter = logging.Formatter('spectre::%(asctime)s::%(levelname)s::%(name)s>  %(message)s')
    if not py_logger.handlers:
        py_handler = logging.StreamHandler(sys.stderr)
        py_handler.setFormatter(py_formatter)
        py_logger.addHandler(py_handler)
    return py_logger
