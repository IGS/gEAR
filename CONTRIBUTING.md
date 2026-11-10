# Contributing to gEAR

> **Draft:** This is a first pass and is under review by the gEAR team.

Thank you for your interest in gEAR. We welcome bug reports, documentation fixes, feedback on the portals and code contributions.

Everyone who takes part in the project is expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

* **Report a bug or ask a question.** See [Reporting bugs and asking questions](#reporting-bugs-and-asking-questions).
* **Suggest a feature or give feedback.** Open a GitHub issue or use the Contact form on any gEAR portal.
* **Improve the documentation.** The developer and analyst docs live in [docs/](docs/README.md), and the user guide lives in [docs/wiki](docs/wiki/gEARWiki.md). Small fixes are a good first contribution.
* **Contribute code.** See [Making a code change](#making-a-code-change).

## Reporting bugs and asking questions

1. Search the [issue tracker](https://github.com/IGS/gEAR/issues) to see whether the problem has already been reported.
2. If not, [open a new issue](https://github.com/IGS/gEAR/issues/new). Please include:
   * the portal URL (for example <https://umgear.org>) and the page you were on
   * the gene, dataset, gene list, etc involved, if any
   * Any parameters used for the analysis or to generate the view
   * what you expected to happen and what happened instead
   * your browser and operating system
   * a screenshot, if it helps

You can also use the Contact form on any gEAR portal. It creates a GitHub issue for you. If your message involves private data, such as an unpublished dataset or your email address, tick the private box so it goes to a private tracker instead of the public one.

**Security problems:** do not open a public issue. Follow the steps in [SECURITY.md](SECURITY.md).

## Making a code change

### Set up a development environment

Running a full gEAR instance takes some work. Start with:

* the [Developer Guide](docs/developer/README.md) for the architecture and workflow
* the [Setup Guides](docs/developer/setup/README.md) for server installation, including a Docker environment
* the [Code Map](docs/developer/code_map.md) to find your way around the code

Then configure the instance and install the git hooks:

```bash
cp gear.ini.template gear.ini   # then edit gear.ini for your setup
pip install pre-commit
pre-commit install
```

The pre-commit hook bumps `www/cache_version.json` on each commit so browsers pick up new CSS and JS. See the [cache busting guide](docs/developer/misc/cache_busting_guide.md).

### Branches and pull requests

We follow the [nvie git branching model](https://nvie.com/posts/a-successful-git-branching-model/):

* `main` is production code
* `devel` is where development happens
* create feature branches from `devel` as `feature/<short-name>`, and hotfix branches as `hotfix/<short-description>`

If you are not on the gEAR team, fork the repository and create your branch in your fork.

When your change is ready, open a pull request against `devel`. Please:

* keep each pull request to one topic
* link the issue it addresses (for example, "Fixes #123")
* describe what you changed and how you tested it
* include screenshots for user interface changes

A maintainer will review the pull request and may ask for changes before merging it.

### Code style

Follow the [code style guidelines](docs/developer/README.md#code-style-guidelines) in the Developer Guide. For Python, run [Ruff](https://docs.astral.sh/ruff/) with the repository's `ruff.toml`:

```bash
ruff check .
```

### Tests

The server-side test suite uses pytest. CI runs it on pull requests to `devel` and `main` that touch server code, services or docs. It does not need a live server or database:

```bash
pip install -r tests/python/requirements.txt
cd tests/python && python -m pytest
```

Please add or update tests for your change. A bug fix should come with a regression test that fails without the fix. CI also runs CodeQL and `pip-audit`. There is no maintained automated test suite for the browser code yet, so please describe how you tested any JavaScript or page changes in your pull request. See [testing.md](docs/developer/testing.md) for details.

### Documentation

If your change affects how gEAR is installed, configured or used, update the matching page under [docs/](docs/README.md). The test suite checks that relative links in the Markdown docs resolve. Do not edit `docs/posters` or `docs/ui-v2-design`, which are kept for historical reference.

## Citing gEAR

If you use gEAR in your research, please cite it as described in [CITATION.cff](CITATION.cff). GitHub also shows the citation under "Cite this repository" on the repository page.

## License

gEAR is released under the [GNU Affero General Public License v3.0](LICENSE). By contributing, you agree that your contributions will be released under the same license.
