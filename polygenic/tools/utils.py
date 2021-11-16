import os
import logging
import sys

def error_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def error_exit(e):
    time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    error_print(" -> polygenic " + version + " failed at " + time)
    error_print(" -> with command: pgstk " + sys.argv.join(" "))
    error_print(" -> with message: ")
    error_print(str(e))
    exit(1)

def expand_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path)) if path else ''

def setup_logger(path):
    logger = logging.getLogger('pgstk')

    log_directory = os.path.dirname(os.path.abspath(os.path.expanduser(path)))
    if log_directory:
        try:
            os.makedirs(log_directory)
        except OSError:
            pass
    logger.setLevel(logging.DEBUG)
    logging_file_handler = logging.FileHandler(path)
    logging_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logging_file_handler.setFormatter(formatter)
    logger.addHandler(logging_file_handler)

    return logger