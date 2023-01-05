===========
Development
===========

Code development follows a collaborative model of forks_ (or clone) and `pull requests`_
on https://github.com using Git_ (see the :ref:`Basic info on git section<basis-info-on-git>`).

.. _Git: https://git-scm.com/
.. _`forks`: https://docs.github.com/en/get-started/quickstart/fork-a-repo
.. _`pull requests`: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests

* Make a new branch for whatever you are doing.  When you are done, push
  it to your own repository and make a merge request from that branch in your
  repository to official master.


Forking
-------

Forking (see forks_) QuantumTransport give you your own personal copy of 
the project at
https://github.com/your-user-name/QuantumTransport

  ``origin`` refers to your copy of the project located at
    https://github.com/your-user-name/QuantumTransport


* Clone your fork to your local machine::

      $ git clone https://github.com/your-user-name/QuantumTransport
      $ cd QuantumTransport
      $ git remote -v

  List the remotes that you are tracking.

Cloning
-------
* Clone QuantumTransport to your local machine and change ``origin``::

      $ git clone https://github.com/nguyenka/QuantumTransport
      $ git remote add orgin https://github.com/your-user-name/QuantumTransport
      $ git remote -v

  List the remotes that you are tracking should have a new ``origin``.


Making changes
--------------

Changes and/or additions can be made both directly in Github for small
changes (see the :ref:`small changes section<making-small-changes>`) and on a
local branch in your fork.  

  * First, checkout a (new) local branch. For example, 
    to make changes to the file develop.rst:

        $ git checkout -b add-develop-rst

  * After making changes, stage the file to be committed using ``git add``::

        $ git add develop.rst

  * Check status of the repos::

        $ git status

  * Commit the staged changes and add commit message or commit text file::

        $ git commit -m "ENH: Add developer workflow guidelines"
        $ git commit -a -F commit.txt
        $ git commit

    A text editor will appear using ``git commit``. Type in commit message
    using guidelines from `here <http://chris.beams.io/posts/git-commit/>`_:

    The seven rules of a great git commit message

      1. Separate subject from body with a blank line
      2. Limit the subject line to 50 characters
      3. Capitalize the subject line
      4. Do not end the subject line with a period
      5. Use the imperative mood in the subject line
      6. Wrap the body at 72 characters
      7. Use the body to explain what and why vs. how

    Read the :ref:`commit message
    section<writing-the-commit-message>` guidelines for commit messages for
    additional information.

  * Push commits to your Github repository::

        $ git push -u origin add-develop-rst

  * `pull requests`_ and wait for feedback and address concerns as needed
    by adding more commits to the 'add-develop-rst' branch on your 
    personal repository and then pushing to your github repository.

  * After the merge-request is done, delete the branch locally::

        $ git branch -D add-develop-rst


.. _making-small-changes:

Making small changes
--------------------

If you want to fix a typo or a line of code. Github has a web_ editing feature. 
Here are the steps to do that there:

* Go to https://github.com/your-user-name/QuantumTransport
* In your repository, browse to the file you want to edit and fix the typo or code
* When you edit a file in another user's repository, github will automatically fork the repository and open a pull request for you.

.. _web: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files

.. _writing-the-commit-message:

Writing the commit message
--------------------------

Here are the ollowing the guidlines from ASE, commit messages should be clear and follow a few basic rules.  Example::

   ENH: add functionality X to QuantumTransport.<submodule>.

   The first line of the commit message starts with a capitalized acronym
   (options listed below) indicating what type of commit this is.  Then a blank
   line, then more text if needed.  Lines shouldn't be longer than 72
   characters.  If the commit is related to a ticket, indicate that with
   "See #3456", "See ticket 3456", "Closes #3456" or similar.

Describing the motivation for a change, the nature of a bug for bug fixes or
some details on what an enhancement does are also good to include in a commit
message.  Messages should be understandable without looking at the code
changes.  A commit message like ``MAINT: fixed another one`` is an example of
what not to do; the reader has to go look for context elsewhere.

Standard acronyms to start the commit message with are:

:API: API change
:BUG: bug fix
:DEP: deprecate something, or remove a deprecated object
:DEV: development tool or utility
:DOC: documentation
:ENH: enhancement
:MAINT: maintenance commit (refactoring, typos, etc.)
:REV: revert an earlier commit
:STY: style fix (whitespace, PEP8)
:TST: addition or modification of tests

Code review
===========

Before you start working on a Merge Request, *please* read
:ref:`coding conventions`.  Please also install pylint!

Hopefully someone will look at your changes and give you some
feedback.  Maybe everything is fine and things can be merged to the official
repository right away, but there could also be some more work to do like:

* make it compatible with Pythons 3.6 and 3.8 
* write more comments
* fix docstrings
* write a test
* add some documentation

.. _basis-info-on-git:

Basic info on git
------------------

* `Git Reference <http://gitref.org>`__
* `Pro Git <https://git-scm.com/book/en/v2>`__
* `Introduction to Git with Scott Chacon of GitHub
  <https://www.youtube.com/watch?v=ZDR433b0HJY>`__
* `Tech Talk: Linus Torvalds on git
  <https://www.youtube.com/watch?v=4XpnKHJAok8>`__

