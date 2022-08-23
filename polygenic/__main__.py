# -*- coding: utf-8 -*-

"""
Entrypoint module
"""

import sys

from polygenic.pgstk import main

if __name__ == "__main__":
    main(sys.argv[1:])