# kaepora

Type Ia Supernova spectral database and tools for creating composite spectra.


## Getting Started



### Prerequisites

Python 2.7

numpy
matplotlib
sqlite3
scipy
astropy
specutils

Version specific dependencies:
msgpack-python version 0.4.6
msgpack-numpy version 0.3.5 or 0.3.6

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

kaepora
=======
Spectral template code using an SQL database for Type Ia Supernova observations.

Version specific dependencies:
msgpack-python version 0.4.6
msgpack-numpy version 0.3.5 or 0.3.6

Example Query:
python query_db.py nb "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN and phase between -1 and 1"

	First argument "b" additionally estimates errors via bootstrap resampling for template spectra. "nb" generates template spectra without errors (much faster).