"""
pdefilepath.py
-----------------
Stores the file path for the NASA JPL ephemeris data used in astronomical computations.

This file defines the path to the ephemeris binary (e.g., DE421, DE441) in a
way that works reliably regardless of the current working directory.
"""

from pathlib import Path

# Define the path to the ephemeris file relative to this script's location
# This ensures the path works even if the script is run from a different folder
EPHEM_PATH: str = str(Path(__file__).parent / "ephem" / "de441.bsp")

# Example usage:
# from pdefilepath import EPHEM_PATH
# planets = load(EPHEM_PATH)
