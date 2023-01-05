.. _coding conventions:

==================
Coding Conventions
==================

Importing modules
=================

Do *not* use the ``import *`` syntax.  Explicitly import everything::

  from munch import Munch
  from cclib.parser import utils


Python Coding Conventions
=========================

Using `pep8 or pyflint` will help improve the readability and consistency of
Python code, as described below.

Use
---
* Use "StudlyCaps" for class names.
* Use 4 spaces for each indentation level.  No hard tabs.
* Use lowercase or "lowercase_with_underscores" for function,
  method, and variable names.  
* Use 'single quotes' for strings, and """triple 
  double quotes""" for ``docstrings``.  Double quotes are OK for
  something like ``"don't"``.

Avoid
-----
* Avoid trailing whitespaces.
* Avoid one-liner compound statements (e.g., ``if x: return``)
* Avoid lambda expressions, which are difficult to understand. 
* Void lines more than 78 characters long.


.. _Style Guide for Python Code:
.. _PEP8: http://www.python.org/peps/pep-0008.html
.. _Docstring Conventions: http://www.python.org/peps/pep-0257.html
.. _Docutils project: http://docutils.sourceforge.net/docs/dev/policies.html
                      #python-coding-conventions
.. _trailing whitespaces: http://www.gnu.org/software/emacs/manual/html_node/
                          emacs/Useless-Whitespace.html


Writing documentation in the code
=================================

Here is an example of how to write good docstrings:

  https://github.com/numpy/numpy/blob/master/doc/example.py

Run pep8 and pylint on your code
================================

It's a good idea to run both the `pep8
<http://pep8.readthedocs.org/en/latest/index.html>`__ and pylint on
your code (or use a text editor that does it automatically)::

    $ pep8 --ignore W293,E129 filename.py
    $ pylint filename.py


