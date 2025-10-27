# An indirect approach to calculate spatial correlation of AvgSa(T)
In this repository you can find some scripts that provide a helpful calculation of an indirect approach of AvgSa(T), and some comparisons when used different inter-Sa(T) spatial correalation models, comparison to a direct formulation of spatial correlation modelling of AvgSa(T). 
Also comparisons using within-event residuals and total-residuals is provided herein.

For more detail please see:
# Reference
Monteiro, V.A, and O’Reilly, G.J. (2026) ‘Notes on spatial correlation for average spectral acceleration: direct and indirect approaches’, (under review)

# How to use
### Define the IM inter-distance. See example.py

```python
---------------  MAO2026 function  ---------------
Arguments:
    IM1 (float): First intensity measure.
    IM2 (float): Second intensity measure.
    h (float): 
    cluster (int, optional): Clustering flag. Defaults to 0.
        - 0 → use "non_cluster" database (ignores vs30).
        - 1 → use "cluster" database (requires vs30).
    vs30 (int, optional): Site condition index. Used only if cluster=1.
        - 1 → low Vs30 (soft soil).
        - 2 → high Vs30 (stiff soil).
        If cluster=0, this is ignored.
Returns:
    corr (float): Correlation value retrieved from the chosen database,
                    for the given IM1, IM2, and distance h.
                    
IMs available are 'Sa', 'Saavg2', 'Saavg3', 'FIV3', 'PGA', and 'PGV'.
For Sa the range of periods is [0.01, 5.0]s
For Saavg2, Saavg3, and FIV3 the range of periods is [0.1, 4.0]s


from GetSpatialCorrelation import MAO2026

# Example
IM1 = "FIV3(1.5)"
IM2 = "Sa(0.1)"
h_distance = 20 #km
corr = MAO2026(IM1, IM2, h_distance, cluster=1, vs30=2)
