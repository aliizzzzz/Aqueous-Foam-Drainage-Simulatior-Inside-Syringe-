# Foam-Drainage-Simulator
This script needs to be run via an IPython console for the most optimum user experience.

## Disclaimer
**This simulation gives the user the option to select two sclerosing foam formulation techniques (the Double Syringe System and the Tessari method) and different liquid-to-gas ratios (1:3, 1:4 and 1:5), however, all formulations use the drainage kinetics of 1:4 DSS foams due to lack of data on other formulations. Data for drainage kinetics of other foam formulations will be added in the near future**

## The Science:
[Sclerosing Foams](https://www.riaendovascular.com/services/foam-sclerotherapy/) are mixtures of a surfactant solution with room air that are produced in-house by clinicians and are administered to treat [varicose veins](https://en.wikipedia.org/wiki/Varicose_veins).

As mixtured of liquid and gas, sclerosing foams are thermodynamically unstable; meaning as time passes by, the liquid content draings and bubbles coarse and coalesce together. One way to quantify foam stability, is by calculating its "half-life" (the time required for half of the liquid content of the foam to drain). In a [previous article](https://journals.sagepub.com/doi/full/10.1177/0268355515589063#_i15), we have characterised the half-life of sclerosing foams manufactured using the Double Syringe System (DSS) with a 1:4 liquid-to-gas ratio. 

Typically, sclerosing foams are injected into the patient's veins using a syringe. An gap of information in the literature is:
  
## Does the act of injection affect the rate of drainage?

The bench-top "drainage" rate of sclerosing foams is characterised previously in the literature (although not extensively). However, their drainage rate during injection, inside the syringe, remains uncharacterised. So the question is: Does the act of injection affect the rate at which the liquid content of a given slcerosing foam drains? So, This script can simulate the injection of foam using a horizontal syringe and outputs a graph of drained liquid height over time. This model uses values of drainage rate for foams that are not under flow. So, if we record vidoes of foam injection, and use image processing to obtain liquid height profiles, we can compare experimental profiles with simulated profiles. If they match, then the act of injection does not contribute to the drainage phenomenon. If they do not match (given a degree of statistical significance), then the act of injection contributes to the rate of drainage and treatment may be improved by employing lower injection rates.

Mathematical equations governing this phenomenon will be added as supplementary material soon!

## Acknowledgements
This is a side-project developed by PhD candidate [Alireze Meghdadi](https://www.linkedin.com/in/alirezameghdadi/).

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.
