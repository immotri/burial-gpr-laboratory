# Contributing

Use an issue to describe a reproducible problem or modelling proposal. Include
the project JSON, expected behaviour, observed behaviour and software versions.
Use synthetic examples that others can reproduce.

For a code change:

1. Keep geometry and material calculations consistent in Python and JavaScript.
2. Run the checks listed in the README.
3. For a numerical change, explain its physical assumptions and validation.
4. Preserve physical units, random seeds, actual survey coordinates and metadata.
5. Keep numerical convergence distinct from a successful solver run.

`toolkit/anatomy.json` and `web/anatomy-data.mjs` are generated together by
`python scripts/build_anatomy.py`. Do not edit just one copy.

A full-wave calculation is not required for documentation-only edits. Changes
to model export should be checked with an appropriately sized gprMax case.
