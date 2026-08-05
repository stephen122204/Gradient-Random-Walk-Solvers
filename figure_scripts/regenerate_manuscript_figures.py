#!/usr/bin/env python3
"""Backward-compatibility wrapper. The canonical script is regenerate_paper_figures.py."""

from regenerate_paper_figures import *  # noqa: F401,F403
from regenerate_paper_figures import main

if __name__ == "__main__":
    main()
