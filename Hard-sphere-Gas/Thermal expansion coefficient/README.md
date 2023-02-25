The thermal expansion coefficient ($\alpha$) can be calculated as:

$$
    \alpha = \frac{1}{V}\left( \frac{\partial V}{\partial T} \right)_P= -\frac{1}{V} \left( \frac{\partial V}{\partial P}\right)_T \left( \frac{\partial P}{\partial T}\right)_V = -\left(\frac{\partial \ln{V}}{\partial P}\right)_T \cdot \left( \frac{\partial P}{\partial T}\right)_V
$$

Where the change of variables is used because in this simulation, the pairs temperature-pressure and pressure-volume are much easier to implement simultaneously as required than volume-temperature.

The file <code> regressions4coefficient.py </code> shows the variation of pressure with respect to temperature and the variaton of the logarithm of the volume with respect to the pressure. It also calculates the slope of the two plots, which are the parameters that, when multiplied, give the numerical thermal expansion coefficient.
