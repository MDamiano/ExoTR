# ExoTR

Version 1.2.0

A fully Bayesian inverse retrieval algorithm to interpret exoplanetary transmission spectra.

Includes:
* A generation of transmission spectra routine;
* A retrieval routine based on nested sampling (i.e. MultiNest).

## Authors
* Mario Damiano (Jet Propulsion Laboratory, California Institute of Technology)
* Renyu Hu (Jet Propulsion Laboratory, California Institute of Technology)

## Collaborators
* Armen Tokadjian (Jet Propulsion Laboratory, California Institute of Technology)
* Audrey DeVault (Massachusetts Institute of Technology)

## Installation:
Download the .zip file from Github, unzip it and place it in a folder at your preference. 
Therefore, make the folder searchable for python in your `.bash_profile` or `.bashrc` depending on your system

`export PYTHONPATH="$PYTHONPATH:/full/path/of/folder/containing/ExoTR:"`

Download the folder "Data" from the following link : [Google Drive](https://drive.google.com/drive/folders/1yXtKIHfsHfCrS9kJef0ycjK1zOvz2vXx?usp=share_link) .
Place the downloaded "Data" folder inside the ExoTR folder.

## Usage
You have to prepare the "retrieval_example.dat" and "forward_example.dat" parameters files before running ExoTR. Refer to the examples provided for guidance.
The full list of possible parameters are listed in the "standard_parameters.dat" file, placed inside the ExoTR package. Do not modify the "standard_parameters.dat" file.

You can generate a transmission spectrum by typing in a python instance or script.py file the following lines:

`import ExoTR`

`spec = ExoTR.CREATE_SPECTRUM()`

`spec.run_forward('forward_example.dat')`

You can run a retrieval by typing in a python instance or script.py file the following lines:

`import ExoTR`

`ret = ExoTR.RETRIEVAL()`

`ret.run_retrieval('retrieval_example.dat')`

To run the retrieval mode you need to have the MultiNest libraries installed in your system as well as `pymultinest (v2.11)`.

`pymultinest` is MPI compatible, therefore, you can run ExoTR to perform the sampling of the retrieval in parallel (you will need to install `mpi4py`):

`mpirun -np 10 python exotr_retrieval.py`

## Plotting the results
The plotting of the retrieval results is automatic and will produce the following graphs:
* Chemistry of the atmosphere versus the atmospheric pressure;
* Mean molecular mass versus the atmospheric pressure;
* The input spectral data and the best fit model calculated by the Bayesian sampling;
* The contribution plot;
* The traces of the fitted free parameters;
* The posterior distribution corner plot.

In case `pymultinest` finds multiple solutions, ExoTR will automatically plot the aforementioned graphs for each of the solutions.

## Acknowledgement
The research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004).

## License
Copyright © 2024, by the California Institute of Technology. ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any commercial use must be negotiated with the Office of Technology Transfer at the California Institute of Technology.

This software may be subject to U.S. export control laws. By accepting this software, the user agrees to comply with all applicable U.S. export laws and regulations. User has the responsibility to obtain export licenses, or other export authority as may be required before exporting such information to foreign countries or providing access to foreign persons.

Licensed under the Apache License, Version 2.0 (the "LICENCE");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## BUGS!!!
For any issues and bugs please send an e-mail at [mario.damiano@jpl.nasa.gov](mario.damiano@jpl.nasa.gov), or submit an issue through the Github system.
