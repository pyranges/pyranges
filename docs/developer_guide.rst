Developer guide
===============

PyRanges is mainly developed by Endre Bakken Stovner and by the Comparative Genomics lab of 
Marco Mariotti. It follows the guidelines for open source software, and external contributors 
are welcome. The code is centralized on github, at https://github.com/pyranges/pyranges

Bugs and feature requests can be reported as github issues. You may also contribute by submitting 
your own code edits, either to deal with issues or to add new functionalities. The code will be 
reviewed by the core development team, which may integrate it in the main repository. Contributions 
are tracked by github and are publicly visible.

Below, we sketch a guide to contribute to PyRanges. It assumes familiarity with python and with the 
terminal, and minimal experience with git/github. Before the actual list of steps follow (Task 
sequence), we go over some essential concepts used in the "continuous integration" system in place 
to maintain and evolve pyranges.



Tests
~~~~~

Tests are an essential part of continuous integration. Briefly, they ensure that code edits do not 
break existing functions. Various layers of tests are implemented in pyranges:

- **unit tests**: quick and compulsory tests about the main pyranges functionalities
- **doctest**: quick and compulsory tests that ensures that the code in the documentation (tutorial and how-to-pages) gives the expected results
- **property based tests**: time-consuming tests that involve the generation of random data to check that the results of pyranges functions match that of other reference bioinformatic tools. These tests are not compulsory: the core development team runs them when the code backbone is edited.

If the code submitted to pyranges does not pass the compulsory tests, it will not be integrated. 
Therefore, we highly recommend developers to run tests before code submissions, as explained 
further below.



Documentation: docstrings
~~~~~~~~~~~~~~~~~~~~~~~~~

Python docstrings are widely used to document the rationale, input arguments, and returned values of 
all functions and methods. The use of a consistent docstring style allows the automatic generation 
of API documentation, as seen in pyranges documentation at `https://pyranges.readthedocs.io/ 
<https://pyranges.readthedocs.io/>`_, built through the Sphynx software.

Pyranges adopts the NumPy/SciPy-style: `https://numpydoc.readthedocs.io/en/latest/format.html 
<https://numpydoc.readthedocs.io/en/latest/format.html>`_. It is important that code contributors 
who edit any function also update their docstrings to reflect how it works; and that all new 
functions contain an appropriate docstring. Follow the link above and inspect existing pyranges 
code to write good docstrings.



Code formatting and linting
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pyranges code follows strict guidelines about its formatting and non-redundancy. This burden is 
not upon the developer: instead, this is achieved by running dedicated software that polishes, 
formats, and "lints" the code before its integration in the main repository. These tools are used:

- **black**: code formatter
- **flake8**: code linting
- **isort**: sort import statements
- **mypy**: type checking (optional, for functions with a defined return type)



Task sequence
~~~~~~~~~~~~~


1. Set up your developer environment
------------------------------------

We recommend creating an environment dedicated to the development of pyranges:

.. code:: bash
	
	conda create -n prdev python pip
	conda activate prdev

Next, install additional packages required for testing and linting code, and for checking the 
automated documentation:

.. code:: bash

	pip install flake8 black isort hypothesis pytest \
	 pytest-watch sphinx sphinx-autoapi sphinxcontrib-napoleon \
	 pyfaidx pyBigWig sphinx_rtd_theme


2. Create and setup your own PyRanges fork
------------------------------------------

The easiest way to do this is through github. Login into the github website if you aren't already, 
then visit `https://github.com/pyranges/pyranges <https://github.com/pyranges/pyranges>`_, and 
click "Fork" on the top right. Fill the form and confirm. In the page of your new fork, find the 
**<> Code** button, and copy the https address. On your computer, create a new folder dedicated 
to the project, then clone your fork inside it:

.. code:: bash

	mkdir pr_debug
	cd pr_debug
	git clone PASTE_GITHUB_HTTPS

Next, cd into your pyranges fork, and install it locally with pip as shown below. By using pip 
option ``-e``, your installation is in "editable" mode: any changes you make to your pyranges code 
will be immediately reflected in your environment. In other words, you won't need to re-run pip 
install every time you change something in the code.

.. code:: bash

	cd pyranges
	pip install -e .


3. Edit the code
----------------

Now, you're ready to edit the code in the pyranges/ folder. 

To run your code to see that it behaves as intended, we recommend using a separate script that 
imports pyranges, making sure you're in the prdev conda environment.


4. Format and lint code
-----------------------

Run these commands from inside your pyranges folder (which has a pyranges subfolder):

.. code:: bash

	black -l 120 pyranges/
	isort --profile black -l 120 tests pyranges
	flake8  --max-line-length=120 --ignore E203,E501,W503 tests pyranges


5. Inspect the Sphynx documentation
-----------------------------------

Your code edits may warrant edits in the pyranges docstrings. In this case, it is compelling to 
locally check that the automatically generated documentation is built appropriately. Inside the 
pyranges folder, run these commands:

.. code:: bash 
	
	cd docs
	make html
	cd -

If the "make" command has no major errors, it will generate the full pyranges documentation in the 
form of html pages, identical to `https://pyranges.readthedocs.io/ <https://pyranges.readthedocs.io/>`_. 
Open the file docs/build/html/index.html with a browser to inspect all the parts that may have 
been affected by your changes, and fix any potential problems. To know more about its inner workings, 
read about the Sphynx system.


6. Run tests
------------

For each of the tests, inspect the output of py.test: warnings are acceptable, but errors must be 
fixed. To run the **unit tests**, run this from inside your pyranges folder:

.. code:: bash

	py.test tests/unit

To run the **doctest**, run this:

.. code:: bash

	py.test tests/tutorial_doctest

To run the non-compulsory **property-based tests**, run:

.. code:: bash
	
	py.test tests/property_based/

If all tests have worked correctly, you are ready to submit your code for integration into the 
main pyranges repository; that is to say, to open a "pull request". Before you can do that, you 
have to update your remote repository, i.e. your pyranges fork at github.

7. Log your changes
----------------------

First, bump the version number in the file pyproject.toml. Then, it's essential to document your changes 
in the CHANGE_LOG.txt file. This log should provide a clear and
concise summary of the modifications, additions, and fixes made in each version of your project. Include
relevant details such as feature enhancements, bug fixes, and any other notable alterations to help
maintain a transparent and informative record of your project's evolution.

8. Push to your remote repository
---------------------------------

Run this command to list all the local files you modified:

.. code:: bash 
	
	git status

You must tell git which of these files have to be synchronized, i.e. "git add" them. You can do this 
by explicitly providing the list of files with:

.. code:: bash

	git add file1 file2 ... fileN

Alternatively to the previous command, if you want to add ALL edited files, you can use:

.. code:: bash

	git add . --dry-run

to check the list of all modified files, then this to actually add them:

.. code:: bash

	git add .

After adding files, you have commit your changes locally with:

.. code:: bash

	git commit -m"Include an informative message here"

Finally, push to your remote repository, i.e. update your online fork at github, with:

.. code:: bash

	git push

You will be requested your github credentials. Note that your online password may not work; in this 
case, google how to set up a github token that you can use.


9. Open a pull request
----------------------

The easiest way to open a pull request is through the github website. Go to **your** 
pyranges fork on github, then find the "Contribute" button (near the **<> Code** button). Click 
it, and select **Open pull request**.

In the newly opened page, carefully check that source and destination are correctly selected. The 
Base repository should be pyranges/pyranges (i.e. the main pyranges repo), and the Head repository 
should be your fork. If you worked on non-master git branches, select them here.

In the comments, write a summary of the introduced changes and their rationale, tagging any related 
github issues (i.e. paste their http address). On the rest of the page, you are presented with a 
list of the code edits. When you're ready, click "Open pull request".


10. Core team only: upload to PyPI
---------------------------------

Every now and then, the core development team considers that a new pyranges version should be 
released. To do so:

- Update the version number in the pyproject.toml file
- Find the "Build and upload to PyPI" workflow in the left menu of the github actions at `https://github.com/pyranges/pyranges/actions/ <https://github.com/pyranges/pyranges/actions/>`_
- Click the "Run workflow" button on the right

Next, check that everything worked correctly, by confirming that a new pyranges installation via 
pip selects the new version.

Finally, the pyranges conda package at Bioconda is updated automatically upon pip upload. Check 
that this is updated correctly.
