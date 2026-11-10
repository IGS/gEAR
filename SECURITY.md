# Security Policy

> **Draft:** This is a first pass and is under review by the gEAR team.

## Supported versions

gEAR is deployed continuously rather than shipped as numbered releases. Security fixes are made only to:

* the `main` branch of this repository, and
* the portals run by the gEAR team at the University of Maryland School of Medicine Institute for Genome Sciences (IGS), such as [umgear.org](https://umgear.org), [nemoanalytics.org](https://nemoanalytics.org), and the other instances listed in the [README](README.md).

Older commits and tags do not get fixes. If you run your own gEAR instance, update from `main` to pick up fixes.

## Reporting a vulnerability

**Please do not report security vulnerabilities in public GitHub issues, pull requests or discussions, or through a portal's public Contact form.**

Report them privately through GitHub instead:

1. Go to the [Security tab](https://github.com/IGS/gEAR/security) of this repository.
2. Click **Report a vulnerability**.
3. Fill in the form. Only the gEAR maintainers can see the report.

Helpful details to include:

* the affected portal (URL) or file, page, API endpoint or CGI script
* the steps to reproduce the problem, with a proof of concept if you have one
* what an attacker could do with it (for example, read another user's private datasets or run commands on the server)
* any ideas you have for a fix

## What to expect

* We aim to acknowledge your report within 5 business days.
* We will look into the problem and keep you updated on our progress. If we confirm it, we will work on a fix and let you know when it has been deployed.
* Please give us a reasonable amount of time to fix the problem before you disclose it publicly. We will agree on a disclosure date with you.
* If you would like, we will credit you in the fix or the advisory.

## Scope

In scope:

* the code in this repository
* the portals run by the gEAR team at IGS

Out of scope:

* vulnerabilities in third-party dependencies, which should be reported to the upstream project (we track these with `pip-audit` and CodeQL in CI and will update affected packages)
* gEAR instances run by other organizations, which should be reported to whoever runs them

When testing, please:

* do not run denial-of-service attacks or high-volume automated scans against the live portals
* do not access, change or delete data that belongs to other users; use your own account and datasets
* stop and report as soon as you find a problem, rather than exploring further
