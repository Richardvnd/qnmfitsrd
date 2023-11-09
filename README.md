# qnmfits
Least-squares fitting of quasinormal modes to ringdown waveforms. This package follows from Eliot's ['qnmfits'](https://github.com/eliotfinch/qnmfits/tree/main) but with some additions to the code. 

## Installation

The package can be installed with `pip install .` - this assumes you are in the same directory as the `pyproject.toml` file. All dependencies for SXS waveform analysis should be installed automatically. Currently, the [`gwsurrogate`](https://pypi.org/project/gwsurrogate/) and [`surfinBH`](https://pypi.org/project/surfinBH/) packages (needed for the analysis of surrogate models) are not automatically installed.

Note that when first importing the package, there may be a short delay.

## Linting

Code can be tidied automatically by running 'invoke tidy' in the command line. Note that 'black', 'flake8', and 'isort' need to be installed in your environment first. 
