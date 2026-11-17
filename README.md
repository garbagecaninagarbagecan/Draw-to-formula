# Park Track Designer — MUF0091 Group Application Task

A Desmos-style graphing tool built for the MUF0091 (Mathematics Unit 1) Group Application Task: designing a continuous running track for a 3 km × 2 km park using five different types of functions, with labelled lake, playground and café.

## How to run

Open `index.html` in any modern browser. No installation, no internet connection and no dependencies are required — everything runs locally in a single page.

## Controls

| Action | Control |
| --- | --- |
| Draw a track section | **Left-click and drag** on the graph |
| Move (pan) the view | **Right-click and drag** |
| Zoom | Mouse scroll, or the +/−/⌂ buttons |
| Drag lake / playground / café | Switch to **Edit features** mode, then left-drag |
| Delete a selected feature | Select it in Edit mode, press Delete |

## How it works

1. Draw a curve with the left mouse button. On release the stroke is resampled and fitted (least squares) against every candidate model: linear/quadratic/cubic polynomials, square root, rectangular hyperbola, logarithmic, exponential and sine.
2. The best-fitting model is snapped to "nice" coefficients (3 decimal places) and converted into a proper equation with a domain restriction, shown in the left panel in mathematical notation.
3. Each new section is automatically anchored to the nearest open endpoint of the existing track, so the vertical shift makes the junction **exactly continuous** — the value carried into the next equation is the previous section's rounded endpoint.
4. The checklist verifies the task requirements: all five required function types (Study Area 1: polynomial, rectangular hyperbola, square root; Study Area 2: exponential/logarithmic, sine), continuity at every junction, the track staying inside the 3 km × 2 km park in the first quadrant, and total track length.

## For your report

- **Copy** — every section has a Copy button producing a Desmos-ready LaTeX string, e.g. `y=0.557\\sqrt{x-0.68}+0.684\\left\\{0.68\\le x\\le1.7\\right\\}`; paste it directly into Desmos to reproduce the graph there.
- **Copy all (Desmos)** — copies every section's equation at once.
- **Download PNG** — exports the labelled graph for your report.
- **Load example** — loads a complete, valid 5-section design you can start from.
- Click a section in the left panel to isolate it (for "graph each section separately" screenshots); click again to show the full track.
