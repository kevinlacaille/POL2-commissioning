# POL2-commissioning
Noise analysis on POL-2 and sensitivity of JCMT/SCUBA-2 while POL-2 mounted on the telescope.

POL-2 is a rotating half-wave plate polarimeter mounted on the James Clerk Maxwell Telescope (JCMT) in front of the SCUBA-2 camera (https://en.wikipedia.org/wiki/Polarimeter). POL-2 is able to indirectly measure magnetic fields from dusty star-forming regions within our own galaxy. See the following paper our team published, providing first-light analysis of POL-2 with scientific results: https://iopscience.iop.org/article/10.3847/1538-4357/aa70a0/pdf

These codes import commissioning data from POL-2/SCUBA-2 observations on nearby star-forming regions, planets (e.g., Uranus), and distant quasars, which all have known polarizations. We looked at these observations and provided the following analysis:
1. calculated a noise across the central circle of radius 3′ for each map;
2. calculated the average exposure time per pixel across the central circle of 3′ radius for each map;
3. used these results to derive an NEFD (Noise Equivalent Flux Density) to transmission relation across all the maps;
4. used the average exposure time to derive an empirical relation between the exposure time in the central map area and the elapsed time of the map;
5. combined these relationships to get the empirical relationship between the elapsed time of an observation and its RMS, given a specific transmission.

These relationships have been incorporated into the SCUBA-2 ITC in Hedwig. A total of 82 850μm observations were analysed, producing a Q and U polarization maps for each observation, which then were used to generate un-IP-corrected and un-de-biased PI maps.

The noise in each image is estimated by two independent methods:
1. Variance map method: Measure the standard deviation within a radius of 3′ from the variance map.
2. Emission map method: Measure the standard deviation within a radius of 3′ from the emission map, while excluding a central circle of 30′′ radius (quasars) or 40′′ radius (Uranus) in order to mask out the point source at the centre of each map.

The sensitivity for blank regions is measured by the Noise Equivalent Flux Density (NEFD). The NEFD can then be empirically related to the transmission at the time of observation, where NEFD = noise x (exposure time)^(1/2)

Both sensitivity methods give reasonable and precise NEFDs for the Stokes Q and U maps. With the first-light data, we are able to predict the elapsed time required to yield a desired noise value for the central 3′ radius of a POL-2 map.
