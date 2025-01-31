How to contribute
=================

We gladly accept contributions to plantiSMASH. To help us keep track of things
and ensure your contributions can be supported by us in the long term, there
are a couple of guidelines that we need contributors to follow.

Getting started
---------------

- Submit an issue
- Clearly describe the issue including ways to reproduce if it is a bug
- Give some description of what the feature should achieve if it a new feature
- Fork the repository on GitHub 
- Make sure your git configuration includes the correct user name and email address
  you can check these by running `git config user.name` and `git config user.email`, respectively
- You can set these by running `git config user.name "Your Name"` and `git config user.email your.name@example.com`

Making changes
--------------

- Create a topic branch based on the `main` branch
- You can create these by `git checkout -b my_contribution main`
- Please don't work on the `main` branch directly
- Make commits in logical units
- Please follow the PEP8 guidelines when writing python code
- Check for unnecessary whitespace with `git diff --check`
- Use proper commit message formatting:

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

- Ensure all changes are backed up by the necessary tests
- Run all tests to ensure nothing else broke accidentally

Submitting changes
------------------

- Push your changes to a topic branch in your fork of the repository
- Submit a pull request to the plantiSMASH team repository
- Update your issue to include a link to your pull request if the automatic linking did not work
- The plantiSMASH team aims to provide initial feedback on pull requests within three business days.

  Thank you for your contribution! 
