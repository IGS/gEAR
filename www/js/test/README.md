# Tests

## How to set up packages

Ensure npm is installed, and make sure you are in the `<gear_root>/www/js` directory

```bash
npm install --save-dev mocha chai playwright
```

In the created package.json file, add within the outermost braces:

```json
"scripts": {
  "test": "mocha"
}
```

## About the testing packages

Mocha (https://mochajs.org/) is a widely-used testing framework that has a minimal setup, allows for organization of tests, and is flexible to let you choose your own libraries to extend its functionality. Mocha has interfaces for both Behavior-driven development (BDD) and Test-driven development (TDD).

Playwright (https://playwright.dev/) is a framework used for automated end-to-end testing. You can use it to interact and test on page selectors. It also provides headless browser testing by default and works with multiple browser types, including mobile ones.
