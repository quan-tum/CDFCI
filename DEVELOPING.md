# Developing Guide

This is a basic list of things to notice for developing CDFCI.

## Interested in improving the code?
Please [contact us](mailto:zhe.wang3@duke.edu) first or open an issue.

The current code is the beta version and we are still making big changes. Please contact us and get the latest status of the code, so that we can work more efficiently.

***

## Git Development Cycle
We use [GitHub Flow](https://gitversion.readthedocs.io/en/latest/git-branching-strategies/githubflow/) since the project is small. The steps are

1. Update `master` to latest upstream code.
2. Create a feature branch `git checkout -b feature/description`.
3. Do the feature/work.
4. Push feature branch to `origin`.
5. Create pull request from `origin/ -> upstream/master`.
6. Review, fix raised comments, merge your PR or even better, get someone else to.

The codes in `master` and `feature` should always compile. For codes that are unable to compile, please create a new branch `debug` or `working`, and merge back to `feature` after completion or debugging.

## C++ Standard
We use C++14 standard.

## Directory
- ```src``` The source code
  - ```main.cpp``` the file with the ```main``` function
  - ```lib``` the third-party libraries used
- ```test``` The test code
  - ```src``` the test source code
  - ```data``` the input of test quantum systems