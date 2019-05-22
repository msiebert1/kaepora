===============
Getting Started
===============

It appears that the time has finally come for you to start your adventure! You will encounter many hardships ahead... That is your fate. Don't feel discouraged, even during the toughest times!

Prerequisites
=============

Python 2.7

numpy

matplotlib

sqlite3

scipy

astropy

specutils

tabulate

Version specific dependencies:

msgpack-python version 0.4.6

msgpack-numpy version 0.3.6

Downloading/Building the Database
=================================
The source code for *kaepora* can be found at https://github.com/msiebert1/kaepora. First clone the repository using:

.. code-block:: console

    git clone https://github.com/msiebert1/kaepora.git

We recommend that you download the most recent version of the database from https://msiebert1.github.io/kaepora/. Unzip the folder and place the '.db' file in the /data folder of your repository. 

Alternatively, you can build the database from source. This process runs several scripts that homogenize the raw spectral data and takes several hours. If you wish to do this navigate to the /kaepora/src folder and execute the following command:

.. code-block:: console

    python build_kaepora.py

Once you have completed one of these steps you should be ready to interact with the database.