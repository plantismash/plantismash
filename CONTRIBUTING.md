How to contribute
=================

We gladly accept contributions to plantiSMASH. To help us keep track of things
and ensure your contributions can be supported by us in the long term, there
are a couple of guidelines that we need contributors to follow.

Getting started
---------------

- Make sure to have an account on [Bitbucket](https://bitbucket.org/account/signup/)
- Submit a ticket for your proposed change
- Clearly describe the issue including ways to reproduce if it is a bug
- Give some description of what the feature should achieve if it a new feature
- Fork the repository on Bitbucket
- Make sure your git configuration includes the correct user name and email address
  you can check these by running `git config user.name` and `git config user.email`, respectively
- You can set these by running `git config user.name "Your Name"` and `git config user.email your.name@example.com`

Making changes
--------------

- Create a topic branch based on the `master` branch
- You can create these by `git checkout -b my_contribution master`
- Please don't work on the `master` branch directly
- Make commits in logical units
- Please follow the PEP8 guidelines when writing python code
- Check for unnecessary whitespace with `git diff --check`
- Make sure you use commit messages in the proper format

```
(#1234) component: Use short imperative description

The first line should be a very brief summary of the patch, starting with the
issue number from the issue tracker, and the main component changed by the
patch. This should be followed by a blank line followed by a paragraph (or
more) explaining what the change is about, possibly with a reference to
relevant literature when implementing a new prediction module. You can finish
up by adding a line containing the words "fixes #1234" or "implements #1234" to
have the issue # link up to the pull request automatically.
```

- Make sure all changes are backed up by the necessary tests
- Run all tests to ensure nothing else broke accidentally

Submitting changes
------------------

- Push your changes to a topic branch in your fork of the repository
- Submit a pull request to the plantiSMASH team repository
- Update your ticket to include a link to your pull request if the automatic linking did not work.
- The antiSMASH team tries to at least provide initial comments on your pull request within three business days.
