#!/usr/bin/env python3
"""
env_check.py

A lightweight helper script for the PRJNA728240 RNA-seq pipeline.
It prints basic information about the Python runtime and system environment.
This ensures the conda environment was correctly set up.
"""

import sys
import platform

print("=== Environment Check ===")
print("Python version:", sys.version)
print("Platform:", platform.system(), platform.release())
print("Python executable:", sys.executable)
