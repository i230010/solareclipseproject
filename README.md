# Solar Eclipse Calculations

This project demonstrates **solar eclipse computations** in Python, including:

- Eclipse search over a specified date range
- Besselian elements calculation
- Gamma value computation
- Maximum eclipse location (central eclipses)
- Start and end times of the eclipse
- Eclipse path plotting

The project uses **Skyfield** for astronomical calculations, along with custom modules for Besselian and eclipse-related computations.

## Requirements

- Python 3.14
- JPL Development Ephemerides (DEXXX) files located in the `ephem/` directory

## How to Run

This project supports execution using **uv**, a fast Python package manager and runner.

To run the main eclipse computation script:

```bash
uv run main.py
```

## Notes

Accuracy may drift significantly when computing eclipses more than ±1000 years from the present.

Results have been compared against: https://ytliu.epizy.com/eclipse/solar_general.html
