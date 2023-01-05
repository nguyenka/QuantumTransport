=====================
Writing documentation
=====================

The documentation is generated using Sphinx_.  The documentation is
stored on Github as text files in the `docs` directory using the
reStructuredText_ markup language.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _Sphinx: http://sphinx.pocoo.org


Installing Docutils and Sphinx
==============================

.. highlight:: bash

If you do::

    $ pip install sphinx_rtd_theme --user

and add ``~/.local/bin`` to you :envvar:`PATH` environment variable, then
you should be ready to go.  You may need the following installed, but they
are not required: scipy, matplotlib, povray, dvipng, pdflatex, bibtex,
AUCTex, fontconfig, convert (ImageMagick).


.. _using_sphinx:

Using Sphinx
============

First, you should take a look at the documentation for Sphinx_ and
reStructuredText_.

Then :command:`cd` to the :file:`docs` directory and build the html-pages::

  $ cd ~/QuantumTransport/docs
  $ make

Create a branch for your work, make your changes to the ``.rst`` files, run
:command:`make` again, check the results and if things
look ok, create a *merge request*::

    $ git checkout -b fixdoc
    $ idle index.rst
    $ make html
    $ git commit -am "fix typos ..."
    $ git push -u origin fixdoc


Extensions to Sphinx
====================

.. highlight:: rest

Extensions to Sphinx:

**:math:**

   This role is for inline LaTeX-style math.  Example:
   ``:math:`\sin(x_n^2)``` gives you :math:`\sin(x_n^2)`.  This role
   is actually the default for ASE's documentation, so you should leave
   out the ``:math:`` part like here: ```\sin(x_n^2)```.

**.. math::**

   Write displayed LaTeX-style math.  Example::

     .. math:: \frac{1}{1+x^2}

   gives you:

   .. math:: \frac{1}{1+x^2}


