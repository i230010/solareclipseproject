# Solar Eclipse Calculations

This project demonstrates **solar eclipse computations** in Python, including:

- Eclipse search over a specified date range
- Besselian elements calculation
- Gamma value computation
- Maximum Eclipse location
- Start and end times of the Eclipse
- Local Circumstances of an Eclipse
- Eclipse path plotting

The project uses **Skyfield** for astronomical calculations, along with other modules and custom modules for Besselian and Eclipse-related computations.

## Requirements

- Python 3.x
- JPL Development Ephemerides (DEXXX) files located in the `ephem/` directory (must change on the `pdefilepath.py` on `EPHEM_PATH`)

## How to Run

This project supports execution using **uv**, a fast Python package manager and runner.

To run the main Eclipse computation script (on the folder):

```bash
uv run main.py
```

## Notes

- Accuracy may drift significantly when computing eclipses more than ~±1000 years from the present even if use the vector calculations of skyfield. Results have been compared against: https://ytliu.epizy.com/eclipse/solar_general.html 
- Used AI for getting the start end times of Eclipse and Maximum Eclipse for Partial Solar Eclipse
- The local circumstances and central path of Eclipse is from: https://celestialprogramming.com/MeeusEclipseExamples/otherExamples/
- The Besselian Elements generator is from (with minor changes): https://celestialprogramming.com/eclipsesGeneratingBesselianElements/generatingBesselianPolynomials.html 
- Cartopy may not Eclipse pahs properly
