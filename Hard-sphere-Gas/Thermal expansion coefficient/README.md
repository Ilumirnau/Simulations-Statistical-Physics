The thermal expansion coefficient can be calculated as:

$$
    \alpha = \frac{1}{V}\left( \frac{\partial V}{\partial T} \right)_P= -\frac{1}{V} \left( \frac{\partial V}{\partial P}\right)_T \left( \frac{\partial P}{\partial T}\right)_V = -\left(\frac{\partial \ln{V}}{\partial P}\right)_T \cdot \left( \frac{\partial P}{\partial T}\right)_V
$$

Where the change of variables is used because in this simulation, the pairs temperature-pressure and pressure-volume are much easier to implement simultaneously as required than volume-temperature.
